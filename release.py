#!/usr/bin/env python3

import re
import subprocess
import sys

def get_latest_tag():
    try:
        subprocess.run(['git', 'fetch', '--tags'], check=True)
        # Get the most recent tag sorted by version
        result = subprocess.run(
            ['git', 'tag', '-l', 'v*', '--sort=-v:refname'], 
            capture_output=True, 
            text=True
        )
        tags = result.stdout.strip().split('\n')
        return tags[0] if tags else 'v0.0.0'
    except Exception as e:
        print(f"Error retrieving tags: {e}")
        sys.exit(1)


def bump_version(current_version, bump_type):
    """
    Bump the version based on the specified type.
    
    :param current_version: Current version string (e.g., 'v1.4.10')
    :param bump_type: Type of version bump ('major', 'minor', 'patch')
    :return: New version string
    """
    # Remove the 'v' prefix for processing
    version_match = re.match(r'v(\d+)\.(\d+)\.(\d+)', current_version)
    if not version_match:
        print(f"Invalid version format: {current_version}")
        sys.exit(1)

    major, minor, patch = map(int, version_match.groups())

    if bump_type == 'major':
        major += 1
        minor = 0
        patch = 0
    elif bump_type == 'minor':
        minor += 1
        patch = 0
    elif bump_type == 'patch':
        patch += 1
    else:
        print(f"Invalid bump type: {bump_type}. Use 'major', 'minor', or 'patch'.")
        sys.exit(1)

    return f'v{major}.{minor}.{patch}'


def push_new_tag(new_tag, message):
    try:
        subprocess.run(['git', 'tag', '-a', new_tag, '-m', message], check=True)

        subprocess.run(['git', 'push', 'origin', new_tag], check=True)

        print(f"Successfully created and pushed new tag: {new_tag}")
        print(f"Tag message: {message}")
    except subprocess.CalledProcessError as e:
        print(f"Error pushing tag: {e}")
        sys.exit(1)


def main():
    if len(sys.argv) != 3 or sys.argv[1] not in ['major', 'minor', 'patch']:
        print("Usage: python3 release.py [major|minor|patch] \"Commit message\"")
        print("Example: python3 release.py minor \"Add new features\"")
        sys.exit(1)

    bump_type = sys.argv[1]
    commit_message = sys.argv[2]

    if not commit_message:
        print("Error: Commit message must be provided")
        sys.exit(1)

    latest_tag = get_latest_tag()
    print(f"Current latest tag: {latest_tag}")

    new_tag = bump_version(latest_tag, bump_type)
    print(f"New tag will be: {new_tag}")

    confirm = input(f"Do you want to create and push the tag {new_tag} with message '{commit_message}'? (y/N): ")
    if confirm.lower() != 'y':
        print("Tag creation cancelled.")
        sys.exit(0)

    push_new_tag(new_tag, commit_message)


if __name__ == '__main__':
    main()