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


def extract_semvar(from_version_str):
    version_numbers = first(
        re.findall(r"^[vV]?(\d+)\.?(\d*)\.?(\d*)$", from_version_str)
    )
    if version_numbers == None:
        return None

    major_ver = int(version_numbers[0])
    minor_ver = int(f"{version_numbers[1]}") if version_numbers[1] != "" else 0
    patch_ver = int(f"{version_numbers[2]}") if version_numbers[2] != "" else 0

    return f"{major_ver}.{minor_ver}.{patch_ver}"


def check_version_status(git_repo_from_meta):

    host = git_repo_from_meta["host"]
    org = git_repo_from_meta["org"]
    repo = git_repo_from_meta["repo"]
    current_ver = git_repo_from_meta["semver"]

    sem_versions = []
    if host == "github":
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

    elif host == "gitlab":
        encoded = quote_plus(f"{org}/{repo}")
        r = requests.get(
            f"https://gitlab.com/api/v4/projects/{encoded}/repository/tags"
        )

        if r.status_code != 200:
            LOGGER.warning(f"Failed to get tags for {org}/{repo}: {r.json()}")
            return None
        response_data = r.json()

    else:
        raise f"{host} is not supported!"

    available_versions = [x["name"] for x in response_data]

    if available_versions == []:
        LOGGER.warning(f"No versions available for {host}/{org}/{repo}")
        return None

    LOGGER.debug(f"Available versions for {host}/{org}/{repo}: {available_versions}")

    sem_versions = [extract_semvar(x["name"]) for x in response_data]
    sem_versions = [v for v in sem_versions if v != None]
    if sem_versions == []:
        LOGGER.warning(f"{host}/{org}/{repo} versions do not conform to semver")
        return None

    newer_vers = [v for v in sem_versions if semver.compare(v, current_ver) > 0]

    if newer_vers == []:
        LOGGER.info(
            f"{host}/{org}/{repo} does not have a new version compared to {current_ver}, first/last {available_versions[0]}/{available_versions[-1]}"
        )
        return None

    latest_version = max(newer_vers, key=cmp_to_key(semver.compare))

    LOGGER.info(
        f"{host}/{org}/{repo} has a new version {latest_version} compared to {current_ver}"
    )

    return {**git_repo_from_meta, "latest_version": latest_version}


def get_new_versions_from_meta_paths():

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

    git_repos_from_meta = []
    for meta_path in module_meta_paths:
        meta_text = Path(meta_path).read_text()
        main_path = f"{os.path.dirname(meta_path)}/main.nf"
        main_text = Path(main_path).read_text()
        repo = first(
            re.findall(
                r"\s+tool_dev_url:\s+\"?https://(github|gitlab).com/([\w-]+)/([\w-]+)\"?",
                meta_text,
            )
        )

        if repo == None:
            LOGGER.warning(f"No repo found in {meta_path}")
            continue

        LOGGER.debug(f"{meta_path} repo: {repo}")

        if "mulled-v2" in main_text:
            LOGGER.warning(f"Mulled container found in {main_path}")
            continue

        version = first(re.findall(rf".*/([\w-]+):v?([\.0-9]+)", main_text))

        if version == None:
            LOGGER.warning(f"No version found in {main_path}")
            continue

        LOGGER.debug(f"{main_path} version: {version}")

        semver_version = extract_semvar(version[1])

        if semver_version == None:
            LOGGER.warning(
                f"Version {version} from {main_path} does not conform to semver"
            )
            continue

        repo_data = {
            "meta.yml": meta_path,
            "main.nf": main_path,
            "host": repo[0],
            "org": repo[1],
            "repo": repo[2],
            "tool": version[0],
            "semver": semver_version,
        }

        LOGGER.debug(f"{meta_path} version data: {repo_data}")

        git_repos_from_meta.append(repo_data)

    with_latest = [check_version_status(r) for r in git_repos_from_meta]
    filtered = [v for v in with_latest if v != None]
    grouped = {}

    for v in filtered:
        key = f"{v['host']}/{v['org']}/{v['repo']}"
        if key not in grouped.keys():
            grouped[key] = {}
            grouped[key]["latest_version"] = v["latest_version"]
            grouped[key]["modules"] = []
            grouped[key]["modules"].append(v)
            continue

        grouped[key]["modules"].append(v)
        grouped[key]["latest_version"] = semver.max_ver(
            grouped[key]["latest_version"], v["latest_version"]
        )

    return grouped


def get_repo_issue_titles(repo):
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
        x["title"] for x in response_data if x["title"].startswith("[UPDATE]")
    ]

    return update_issues


def register_issue(host_repo_name, issue_title, tool_data):

    first_module = tool_data["modules"][0]
    host, org, repo = (first_module["host"], first_module["org"], first_module["repo"])

    module_paths = "\n".join(
        [f"- {module['main.nf']}" for module in tool_data["modules"]]
    )

    issue_data = [
        "gh api",
        "--method POST",
        '-H "Accept: application/vnd.github+json"',
        '-H "X-GitHub-Api-Version: 2022-11-28"',
        f"/repos/{host_repo_name}/issues",
        f'-f "title={issue_title}"',
        f'-f "body=- https://{host}.com/{org}/{repo}\n{module_paths}"',
        '-f "labels[]=enhancement"',
    ]

    issue_post = subprocess.run(
        " ".join(issue_data),
        shell=True,
        capture_output=True,
        text=True,
    )

    if issue_post.returncode == 0:
        LOGGER.info(f"Submitted issue: {issue_title}")
    else:
        LOGGER.warning(f"Failed to submitted issue: {issue_title}: {issue_post.stderr}")


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

    git_repos_from_meta = get_new_versions_from_meta_paths()
    update_issue_titles = get_repo_issue_titles(PIPELINE_REPO)

    for key, tool_data in git_repos_from_meta.items():

        new_issue_title = f"[UPDATE] {key.upper()} -> {tool_data['latest_version']}"

        if new_issue_title in update_issue_titles:
            LOGGER.info(f"{new_issue_title} has already been raised!")
            continue

        if args.dry_run:
            LOGGER.info(f"Dry run new issue: {new_issue_title}")
            continue

        register_issue(PIPELINE_REPO, new_issue_title, tool_data)
