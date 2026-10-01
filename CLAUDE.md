## Commit attribution policy

Only humans may be represented as authors or contributors.

Do not add AI tools or models as authors, co-authors, contributors, or
signatories in commit messages or pull-request descriptions. This includes
Claude, Codex, ChatGPT, Gemini, Copilot, Cursor, Aider, and Devin.

Never add `Co-authored-by`, `Generated-by`, `Assisted-by`, or similar
attribution referring to an AI system.

Before creating commits, ensure the repository hooks are installed:

    pre-commit install

Never bypass repository hooks with `--no-verify`.

If the hook rejects a commit, remove the prohibited attribution instead of
bypassing the hook.
