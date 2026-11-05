#!/usr/bin/env python3
"""Reject AI-system attribution in Git commit messages.

This is intended for use from pre-commit's commit-msg hook, which passes the
path to Git's temporary commit-message file as its single argument.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


# Add identity terms here as the policy evolves.  The word boundaries prevent
# an identity such as "devin" from matching an unrelated longer word.
AI_IDENTITIES = (
    "claude",
    "anthropic",
    "chatgpt",
    "openai",
    "codex",
    "copilot",
    "gemini",
    "cursor",
    "aider",
    "devin",
    "noreply@anthropic.com",
    "noreply@openai.com",
)

IDENTITY_PATTERN = re.compile(
    r"(?<![\w@.-])(?:" + "|".join(map(re.escape, AI_IDENTITIES)) + r")(?![\w@.-])",
    re.IGNORECASE,
)

# Git trailers and a small set of plain-language forms that are clearly
# attribution.  Ordinary discussion of AI tools elsewhere in a message is
# deliberately outside this match.
ATTRIBUTION_PATTERN = re.compile(
    r"^\s*(?:"
    r"co-authored-by|contributed-by|reviewed-by|acked-by|signed-off-by|"
    r"generated-by|assisted-by|ai-generated-by|ai-assisted-by"
    r")\s*:|^\s*(?:generated|assisted)\s+(?:with|by)\b",
    re.IGNORECASE,
)


def find_prohibited_lines(message: str) -> list[str]:
    """Return attribution lines that identify a prohibited AI system."""
    return [
        line
        for line in message.splitlines()
        if ATTRIBUTION_PATTERN.search(line) and IDENTITY_PATTERN.search(line)
    ]


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(f"Usage: {Path(argv[0]).name} COMMIT_MESSAGE_FILE", file=sys.stderr)
        return 2

    message_file = Path(argv[1])
    try:
        message = message_file.read_text(encoding="utf-8")
    except OSError as error:
        print(f"Cannot read commit message file '{message_file}': {error}", file=sys.stderr)
        return 2
    except UnicodeDecodeError:
        print(f"Commit message file '{message_file}' is not valid UTF-8.", file=sys.stderr)
        return 2

    prohibited_lines = find_prohibited_lines(message)
    if not prohibited_lines:
        return 0

    print("AI attribution is not permitted in commit messages.", file=sys.stderr)
    print("Offending line(s):", file=sys.stderr)
    for line in prohibited_lines:
        print(f"  {line}", file=sys.stderr)
    print("Remove the AI attribution and commit again.", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
