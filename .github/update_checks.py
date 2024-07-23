#!/usr/bin/env python3

import colorlog
from urllib.parse import quote_plus
from functools import cmp_to_key
from pathlib import Path

import subprocess
import argparse
import requests
import logging
import semver
import json
import re
import os

PIPELINE_REPO = "plant-food-research-open/assemblyqc"


def get_logger():
    formatter = colorlog.ColoredFormatter(
        "%(log_color)s%(levelname)-8s%(reset)s %(blue)s%(message)s",
        datefmt=None,
        reset=True,
        log_colors={
            "DEBUG": "cyan",
            "INFO": "green",
            "WARNING": "yellow",
            "ERROR": "red",
            "CRITICAL": "red,bg_white",
        },
        secondary_log_colors={},
        style="%",
    )

    handler = colorlog.StreamHandler()
    handler.setFormatter(formatter)

    logger = colorlog.getLogger("")
    logger.addHandler(handler)
    return logger


LOGGER = get_logger()


def first(from_list):
    if from_list == []:
        return None

    return from_list[0]


def flatten(xss):
    return [x for xs in xss for x in xs]


def extract_semvar(from_version_str):
    version_numbers = first(
        re.findall(r"^[vV]?(\d+)\.?(\d*)\.?(\d*)$", from_version_str)
    )
    if version_numbers == None:
        return None

    minor_ver = f".{version_numbers[1]}" if version_numbers[1] != "" else ".0"
    patch_ver = f".{version_numbers[2]}" if version_numbers[2] != "" else ".0"

    return f"{version_numbers[0]}{minor_ver}{patch_ver}"


def check_version_status(repo_tuple):
    git, org, repo, tag = repo_tuple

    query_ver = extract_semvar(tag)

    if query_ver == None:
        LOGGER.warning(f"{repo_tuple} does not conform to semver")
        return None

    sem_versions = []
    if git == "github":
        r = subprocess.run(
            " ".join(
                [
                    "gh api",
                    '-H "Accept: application/vnd.github+json"',
                    '-H "X-GitHub-Api-Version: 2022-11-28"',
                    f"/repos/{org}/{repo}/tags",
                ]
            ),
            shell=True,
            capture_output=True,
            text=True,
        )

        if r.returncode != 0:
            LOGGER.warning(f"Failed to get tags for {org}/{repo}: {r.stderr}")
            return None
        response_data = json.loads(r.stdout)

    elif git == "gitlab":
        encoded = quote_plus(f"{org}/{repo}")
        r = requests.get(
            f"https://gitlab.com/api/v4/projects/{encoded}/repository/tags"
        )

        if r.status_code != 200:
            LOGGER.warning(f"Failed to get tags for {org}/{repo}: {r.json()}")
            return None
        response_data = r.json()

    else:
        raise f"{git} is not supported!"

    available_versions = [x["name"] for x in response_data]
    LOGGER.debug(f"Available versions for {repo_tuple}: {available_versions}")

    sem_versions = [extract_semvar(x["name"]) for x in response_data]
    sem_versions = [v for v in sem_versions if v != None]
    if sem_versions == []:
        LOGGER.warning(f"{repo_tuple} versions do not conform to semver")
        return None

    newer_vers = [v for v in sem_versions if semver.compare(v, query_ver) > 0]

    if newer_vers == []:
        LOGGER.debug(f"{repo_tuple} does not have a new version")
        return None

    return (git, org, repo, tag, max(newer_vers, key=cmp_to_key(semver.compare)))


def get_new_versions_from_meta_paths() -> list[tuple[str]]:

    module_meta_paths = (
        subprocess.run(
            "find ./modules -name meta.yml",
            shell=True,
            capture_output=True,
            text=True,
        )
        .stdout.strip()
        .split("\n")
    )

    git_repos_by_meta = []
    for meta_path in module_meta_paths:
        meta_text = Path(meta_path).read_text()
        main_text = Path(f"{os.path.dirname(meta_path)}/main.nf").read_text()
        repos = re.findall(
            r"\s+tool_dev_url: \"?https://(github|gitlab).com/([\w-]+)/([\w-]+)\"?",
            meta_text,
        )
        versions = [
            first(re.findall(rf".*/{repo[2]}:([\.0-9]+)", main_text)) for repo in repos
        ]

        repo_ver = [
            (repo[0], repo[1], repo[2], version)
            for (repo, version) in zip(repos, versions)
            if version != None
        ]

        if repo_ver == []:
            continue

        git_repos_by_meta.append(repo_ver)

    git_repos_by_meta = sorted(list(set(flatten(git_repos_by_meta))))

    new_versions = [check_version_status(r) for r in git_repos_by_meta]

    return sorted(
        list(set([(v[0], v[1], v[2], v[4]) for v in new_versions if v != None]))
    )


def get_repo_issues(repo):
    r = subprocess.run(
        " ".join(
            [
                "gh api",
                '-H "Accept: application/vnd.github+json"',
                '-H "X-GitHub-Api-Version: 2022-11-28"',
                f"/repos/{repo}/issues",
            ]
        ),
        shell=True,
        capture_output=True,
        text=True,
    )

    if r.returncode != 0:
        LOGGER.error(f"Failed to get issues for {repo}: {r.stderr}")
        exit(1)

    response_data = json.loads(r.stdout)

    update_issues = [
        first(re.findall(r"\[UPDATE\] (\w+)/([\w-]+)/([\w-]+) -> (.*)", x["title"]))
        for x in response_data
        if x["title"].startswith("[UPDATE]")
    ]

    return update_issues


def register_issue(host_repo_name, tool_data: list[str]):

    issue_data = [
        "gh api",
        "--method POST",
        '-H "Accept: application/vnd.github+json"',
        '-H "X-GitHub-Api-Version: 2022-11-28"',
        f"/repos/{host_repo_name}/issues",
        f'-f "title=[UPDATE] {tool_data[0]}/{tool_data[1]}/{tool_data[2]} -> {tool_data[3]}"',
        f'-f "body=- https://{tool_data[0].lower()}.com/{tool_data[1].lower()}/{tool_data[2].lower()}"',
        '-f "labels[]=enhancement"',
    ]

    issue_post = subprocess.run(
        " ".join(issue_data),
        shell=True,
        capture_output=True,
        text=True,
    )

    if issue_post.returncode == 0:
        LOGGER.info(f"Submitted issue for: {tool_data}")
    else:
        LOGGER.warning(f"Failed to submitted issue for: {tool_data}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose mode"
    )
    parser.add_argument(
        "-d", "--dry-run", action="store_true", help="Enable dry-run mode"
    )
    args = parser.parse_args()

    if args.verbose:
        LOGGER.setLevel(level=logging.DEBUG)
    else:
        LOGGER.setLevel(level=logging.INFO)

    new_versions = get_new_versions_from_meta_paths()
    repo_issues = get_repo_issues(PIPELINE_REPO)

    new_issues = set(
        [(v[0].upper(), v[1].upper(), v[2].upper(), v[3].upper()) for v in new_versions]
    ) - set(repo_issues)

    new_issues = sorted(list(new_issues))

    if new_issues == []:
        LOGGER.info("No updates are available!")
        exit(0)

    LOGGER.info(f"Going to submit following updates: {new_issues}")

    if args.dry_run:
        LOGGER.info(f"Exiting with dry-run")
        exit(0)

    for issue in new_issues:
        register_issue(PIPELINE_REPO, issue)
