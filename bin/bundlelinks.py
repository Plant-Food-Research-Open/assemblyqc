#!/usr/bin/env python3

import argparse


def log(verbose, message):
    if verbose:
        print(message)


def dir_sign(end, start):
    if abs(int(end) - int(start)) == 0:
        return 1
    sign = (int(end) - int(start)) / abs(int(end) - int(start))
    return sign


def is_same_direction(link, bundle, verbose):
    _, ref_start, ref_end, _, target_start, target_end = link
    (
        _,
        bundle_ref_start,
        bundle_ref_end,
        _,
        bundle_target_start,
        bundle_target_end,
    ) = bundle

    link_ref_dir = dir_sign(ref_end, ref_start)
    link_target_dir = dir_sign(target_end, target_start)

    bundle_ref_dir = dir_sign(bundle_ref_end, bundle_ref_start)
    bundle_target_dir = dir_sign(bundle_target_end, bundle_target_start)

    log(
        verbose,
        f"Compared directions l: {link_ref_dir},{link_target_dir} and b: {bundle_ref_dir},{bundle_target_dir}",
    )

    if link_ref_dir == bundle_ref_dir and link_target_dir == bundle_target_dir:
        return True

    return False


def within_max_gap(link, bundle, max_gap):
    _, ref_start, ref_end, _, target_start, target_end = link
    (
        _,
        bundle_ref_start,
        bundle_ref_end,
        _,
        bundle_target_start,
        bundle_target_end,
    ) = bundle

    ref_start = int(ref_start)
    ref_end = int(ref_end)
    target_start = int(target_start)
    target_end = int(target_end)

    bundle_ref_start = int(bundle_ref_start)
    bundle_ref_end = int(bundle_ref_end)
    bundle_target_start = int(bundle_target_start)
    bundle_target_end = int(bundle_target_end)

    ref_within_max_gap = (
        abs(ref_start - bundle_ref_start) <= max_gap
        or abs(ref_start - bundle_ref_end) <= max_gap
        or abs(ref_end - bundle_ref_start) <= max_gap
        or abs(ref_end - bundle_ref_end) <= max_gap
    )

    target_within_max_gap = (
        abs(target_start - bundle_target_start) <= max_gap
        or abs(target_start - bundle_target_end) <= max_gap
        or abs(target_end - bundle_target_start) <= max_gap
        or abs(target_end - bundle_target_end) <= max_gap
    )

    return ref_within_max_gap and target_within_max_gap


def get_bundle_directions(bundle):
    (
        _,
        bundle_ref_start,
        bundle_ref_end,
        _,
        bundle_target_start,
        bundle_target_end,
    ) = bundle

    return (
        dir_sign(bundle_ref_end, bundle_ref_start),
        dir_sign(bundle_target_end, bundle_target_start),
    )


def add_link_to_bundle(link, bundle, verbose):
    ref_dir, target_dir = get_bundle_directions(bundle)

    log(verbose, f"Bundle directions l: {ref_dir},{target_dir}")

    _, ref_start, ref_end, _, target_start, target_end = link

    ref_start = int(ref_start)
    ref_end = int(ref_end)
    target_start = int(target_start)
    target_end = int(target_end)

    (
        ref,
        bundle_ref_start,
        bundle_ref_end,
        target,
        bundle_target_start,
        bundle_target_end,
    ) = bundle

    bundle_ref_start = int(bundle_ref_start)
    bundle_ref_end = int(bundle_ref_end)
    bundle_target_start = int(bundle_target_start)
    bundle_target_end = int(bundle_target_end)

    updated_bundle_ref_start = None
    updated_bundle_ref_end = None
    updated_bundle_target_start = None
    updated_bundle_target_end = None

    if ref_dir > 0:
        updated_bundle_ref_start = min(ref_start, bundle_ref_start)
        updated_bundle_ref_end = max(ref_end, bundle_ref_end)
    else:
        updated_bundle_ref_start = max(ref_start, bundle_ref_start)
        updated_bundle_ref_end = min(ref_end, bundle_ref_end)

    if target_dir > 0:
        updated_bundle_target_start = min(target_start, bundle_target_start)
        updated_bundle_target_end = max(target_end, bundle_target_end)
    else:
        updated_bundle_target_start = max(target_start, bundle_target_start)
        updated_bundle_target_end = min(target_end, bundle_target_end)

    return [
        ref,
        updated_bundle_ref_start,
        updated_bundle_ref_end,
        target,
        updated_bundle_target_start,
        updated_bundle_target_end,
    ]


def bundle_len(bundle):
    (
        _,
        bundle_ref_start,
        bundle_ref_end,
        _,
        bundle_target_start,
        bundle_target_end,
    ) = bundle

    return min(
        abs(int(bundle_ref_end) - int(bundle_ref_start)),
        abs(int(bundle_target_end) - int(bundle_target_start)),
    )


def bundle_links(input_file, output_file, max_gap, min_bundle_size, verbose):
    bundles = {}
    nlinks = {}
    current_bundle_num = 0

    with open(input_file, "r") as f:
        for line in f:
            link = line.strip().split("\t")
            ref, _, _, target, _, _ = link

            log(verbose, f"Link: {link}")

            link_key = f"{ref}:{target}"
            possible_bundles_to_add_to = [
                k for k in bundles.keys() if k.startswith(link_key)
            ]

            log(verbose, f"Possible bundles are: {possible_bundles_to_add_to}")

            # Check if the current link can be added to an existing bundle
            link_added = False
            for k in possible_bundles_to_add_to:
                bundle = bundles[k]

                log(verbose, f"Checking against bundle: {bundle}")

                if not is_same_direction(link, bundle, verbose):
                    log(verbose, "Bundle and link have different directions")
                    continue

                if not within_max_gap(link, bundle, max_gap):
                    log(verbose, f"Bundle and link are not within {max_gap}")
                    continue

                # Add link to bundle
                updated_bundle = add_link_to_bundle(link, bundle, verbose)
                bundles[k] = updated_bundle
                link_added = True
                nlinks[k] += 1
                log(verbose, f"Added link to bundle")
                log(verbose, f"Updated bundle: {updated_bundle}")
                break

            # Create a new bundle for the current link
            if not link_added:
                current_bundle_num += 1
                bundles[link_key + f":{current_bundle_num}"] = link
                nlinks[link_key + f":{current_bundle_num}"] = 0
                log(verbose, f"Created bundle: {link}")

    # Filter out bundles smaller than the minimum bundle size
    filtered_bundles = {
        key: bundle
        for key, bundle in bundles.items()
        if bundle_len(bundle) >= min_bundle_size
    }

    # Write the filtered bundles to the output file
    with open(output_file, "w") as f:
        for k, bundle in filtered_bundles.items():
            bundle_nlinks = nlinks[k]
            f.write(
                "\t".join([str(v) for v in bundle] + [f"nlinks={bundle_nlinks}"]) + "\n"
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Bundle links",
        epilog="Author: ChatGPT, Usman Rashid",
    )
    parser.add_argument("input_file", help="Input TSV file containing links")
    parser.add_argument("output_file", help="Output TSV file to write bundled links")
    parser.add_argument(
        "--max_gap",
        type=int,
        default=1_000_000,
        help="Maximum gap allowed between links for bundling",
    )
    parser.add_argument(
        "--min_bundle_size",
        type=int,
        default=1_000,
        help="Minimum size of a bundle to retain",
    )

    parser.add_argument("--verbose", help="Print info messages", action="store_true")

    args = parser.parse_args()

    bundle_links(
        args.input_file,
        args.output_file,
        args.max_gap,
        args.min_bundle_size,
        args.verbose,
    )
