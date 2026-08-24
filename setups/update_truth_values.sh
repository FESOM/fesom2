#!/bin/bash
# Update the CI truth values that guard FESOM2 results.
#
# Two modes:
#
#   1. Whole pull request (needs a GitHub token):
#        bash update_truth_values.sh 984
#        bash update_truth_values.sh https://github.com/FESOM/fesom2/pull/984
#
#      Finds every failing workflow on the PR head, downloads the job logs,
#      and rewrites every truth value that moved: the fcheck: blocks in
#      setups/<name>/setup.yml and the inline reference dict in
#      .github/workflows/fesom2_xios.yml.
#
#      Add --dry-run to print the changes without writing them.
#
#   Getting a token that works
#   --------------------------
#   Job logs cannot be downloaded anonymously (HTTP 403), so this mode needs a
#   token. Despite GitHub answering "Must have admin rights to Repository",
#   admin is NOT required - read access to Actions is enough.
#
#   github.com -> Settings -> Developer settings -> Personal access tokens
#
#     fine-grained  Resource owner MUST be FESOM, not your personal account -
#                   a token owned by your account cannot read FESOM's logs.
#                   Repository access: only FESOM/fesom2.
#                   Permissions: Actions = Read-only (Metadata is added for you).
#                   Because the owner is an organisation the token may sit in
#                   "pending" until a FESOM owner approves it; until then every
#                   log read 403s even though the scope is correct.
#
#     classic       Scope: repo (public_repo is enough for a public repo).
#                   Blunter, but skips the org-approval step.
#
#     in a workflow secrets.GITHUB_TOKEN already carries actions: read, so this
#                   script can run in CI with no PAT at all.
#
#   Storing it. Keep the token out of the repository - nothing here is
#   gitignored. A private file outside the tree works well:
#
#     printf 'GITHUB_TOKEN=github_pat_...\n' > ~/.ssh/fesom_truth_update_token
#     chmod 600 ~/.ssh/fesom_truth_update_token
#     source ~/.ssh/fesom_truth_update_token
#
#   Write it as a shell assignment as above, not as a bare token. Sourcing a
#   file that contains only the token makes bash try to execute it, and the
#   resulting "github_pat_...: command not found" prints the secret into your
#   terminal and shell history.
#
#   If every log read fails with 401 "Server failed to authenticate the
#   request", that is not the token: see _StripAuthOnRedirect below.
#
#   2. Single setup from pasted CI output (original behaviour, unchanged):
#        bash update_truth_values.sh ci_output.txt [setup.yml]
#        bash update_truth_values.sh setup.yml ci_output.txt
#        bash update_truth_values.sh < ci_output.txt
#        bash update_truth_values.sh setup.yml "paste CI output here"

set -o pipefail

case "${1:-}" in
    -h|--help|help)
        sed -n '2,58p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
        exit 0
        ;;
esac

# ---------------------------------------------------------------- PR mode ---
# Anything that looks like a bare PR number or a pull request URL goes here.
PR_NUMBER=""
DRY_RUN=0
for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=1 ;;
        *[!0-9]*)
            if [[ "$arg" =~ ^https?://github\.com/[^/]+/[^/]+/pull/([0-9]+) ]]; then
                PR_NUMBER="${BASH_REMATCH[1]}"
            fi
            ;;
        ?*) PR_NUMBER="$arg" ;;
    esac
done

if [ -n "$PR_NUMBER" ]; then
    REPO_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
    TOKEN="${GITHUB_TOKEN:-${GH_TOKEN:-}}"
    if [ -z "$TOKEN" ]; then
        echo "Error: job logs cannot be downloaded anonymously, so PR mode needs"
        echo "       a token. See 'Getting a token that works' at the top of this"
        echo "       script, then: source ~/.ssh/fesom_truth_update_token"
        echo
        echo "       No token? The original paste mode needs none - open the failing"
        echo "       job on github.com, copy its output, and run:"
        echo "         $0 setups/<name>/setup.yml \"<paste>\""
        exit 1
    fi

    SLUG=$(git -C "$REPO_ROOT" remote get-url origin 2>/dev/null \
           | sed -E 's|.*github\.com[:/]||; s|\.git$||')
    SLUG="${SLUG:-FESOM/fesom2}"

    DRY_RUN=$DRY_RUN REPO_ROOT="$REPO_ROOT" SLUG="$SLUG" TOKEN="$TOKEN" \
    PR_NUMBER="$PR_NUMBER" python3 - <<'PY'
import json, os, re, sys, urllib.request, gzip, io

root   = os.environ["REPO_ROOT"]
slug   = os.environ["SLUG"]
token  = os.environ["TOKEN"]
pr     = os.environ["PR_NUMBER"]
dry    = os.environ["DRY_RUN"] == "1"

class _StripAuthOnRedirect(urllib.request.HTTPRedirectHandler):
    """Make the log download work.

    /actions/jobs/<id>/logs does not return the log: it 302s to a short-lived
    blob storage URL. urllib follows the redirect but re-sends the GitHub
    Authorization header to that host, which rejects it with

        HTTP Error 401: Server failed to authenticate the request.
        Please refer to the information in the www-authenticate header.

    The message points at the token, but the token is fine - the blob URL is
    pre-signed and must be fetched with no Authorization header at all. So
    drop the header when following a redirect."""
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        new = super().redirect_request(req, fp, code, msg, headers, newurl)
        if new is not None:
            new.remove_header("Authorization")
        return new

_opener = urllib.request.build_opener(_StripAuthOnRedirect)

def api(url, raw=False):
    req = urllib.request.Request(url, headers={
        "Accept": "application/vnd.github+json",
        "Authorization": f"Bearer {token}",
        "User-Agent": "fesom-update-truth",
    })
    with _opener.open(req) as r:
        data = r.read()
    if raw:
        if data[:2] == b"\x1f\x8b":
            data = gzip.decompress(data)
        return data.decode("utf-8", "replace")
    return json.loads(data)

# --- which setup.yml does each workflow drive? ask the workflow itself ------
wf_dir = os.path.join(root, ".github", "workflows")
setup_of = {}
for fn in sorted(os.listdir(wf_dir)):
    if not fn.endswith((".yml", ".yaml")):
        continue
    text = open(os.path.join(wf_dir, fn)).read()
    m = re.search(r"mkrun\s+\S+\s+(\S+)", text)
    if m:
        setup_of[fn] = m.group(1)

print(f"repo   : {slug}")
print(f"PR     : #{pr}")

head = api(f"https://api.github.com/repos/{slug}/pulls/{pr}")["head"]["sha"]
print(f"head   : {head[:10]}\n")

runs = api(f"https://api.github.com/repos/{slug}/actions/runs"
           f"?head_sha={head}&per_page=100")["workflow_runs"]
failed = [r for r in runs if r["conclusion"] == "failure"]
if not failed:
    print("No failing workflows on this head - nothing to update.")
    sys.exit(0)
print(f"failing workflows: {len(failed)}")

# name -> {var: value}, collected from the logs
fcheck_updates = {}   # workflow file -> {var: val}
xios_updates   = {}   # var -> val

RE_FCHECK = re.compile(r"Variable:\s*([^,]+),\s*current_value:\s*([^\s,]+)")
RE_XIOS   = re.compile(r"FAIL:\s*(\w+)\s+got=([-\d.eE+]+)")

for run in failed:
    wf_file = os.path.basename(run.get("path", ""))
    try:
        jobs = api(f"https://api.github.com/repos/{slug}/actions/runs/{run['id']}/jobs")["jobs"]
    except Exception as e:
        print(f"  ! {wf_file}: cannot list jobs ({e})")
        continue
    for job in jobs:
        if job["conclusion"] != "failure":
            continue
        try:
            log = api(f"https://api.github.com/repos/{slug}/actions/jobs/{job['id']}/logs", raw=True)
        except Exception as e:
            print(f"  ! {wf_file}/{job['name']}: cannot read log ({e})")
            continue

        got_f = {v.strip(): n.strip() for v, n in RE_FCHECK.findall(log)}
        got_x = dict(RE_XIOS.findall(log))

        if got_x:
            xios_updates.update(got_x)
            print(f"  {wf_file}: {len(got_x)} XIOS reference value(s) moved")
        if got_f:
            fcheck_updates.setdefault(wf_file, {}).update(got_f)
            print(f"  {wf_file}: {len(got_f)} fcheck value(s) moved")
        if not got_f and not got_x:
            print(f"  {wf_file}/{job['name']}: failed, but no truth mismatch found "
                  f"(different failure - left alone)")

changed = []

# --- rewrite fcheck: blocks in setups/<name>/setup.yml ----------------------
for wf_file, vals in fcheck_updates.items():
    setup = setup_of.get(wf_file)
    if not setup:
        print(f"  ! {wf_file}: no 'mkrun <mesh> <setup>' found, cannot map to a setup.yml")
        continue
    path = os.path.join(root, "setups", setup, "setup.yml")
    if not os.path.exists(path):
        print(f"  ! {path} does not exist")
        continue
    lines = open(path).read().split("\n")
    # only touch keys inside the fcheck: block, so short names like u/v
    # cannot collide with unrelated keys elsewhere in the file
    try:
        start = next(i for i, l in enumerate(lines) if l.strip() == "fcheck:")
    except StopIteration:
        print(f"  ! {path} has no fcheck: block")
        continue
    n = 0
    for i in range(start + 1, len(lines)):
        stripped = lines[i].strip()
        if stripped and not lines[i][:1].isspace():
            break
        m = re.match(r"^(\s*)([A-Za-z_]\w*):\s*(.*)$", lines[i])
        if m and m.group(2) in vals:
            new = vals[m.group(2)]
            if m.group(3).strip() != new:
                lines[i] = f"{m.group(1)}{m.group(2)}: {new}"
                n += 1
    if n:
        changed.append((os.path.relpath(path, root), n))
        if not dry:
            open(path, "w").write("\n".join(lines))

# --- rewrite the inline reference dict in fesom2_xios.yml ------------------
if xios_updates:
    path = os.path.join(wf_dir, "fesom2_xios.yml")
    lines = open(path).read().split("\n")
    n = 0
    for i, l in enumerate(lines):
        m = re.match(r'^(\s*)"(\w+)":(\s*)([-\d.eE+]+),\s*$', l)
        if m and m.group(2) in xios_updates:
            new = xios_updates[m.group(2)]
            if m.group(4) != new:
                lines[i] = f'{m.group(1)}"{m.group(2)}":{m.group(3)}{new},'
                n += 1
    if n:
        changed.append((os.path.relpath(path, root), n))
        if not dry:
            open(path, "w").write("\n".join(lines))

print()
if not changed:
    print("Nothing to update - the failures carry no truth mismatches.")
    sys.exit(0)
for f, n in changed:
    print(f"  {'would update' if dry else 'updated'} {f}  ({n} value(s))")
print()
print("Review with: git diff")
PY
    exit $?
fi

# ------------------------------------------------- single-setup mode -------
# (unchanged from the original script)

# Smart argument detection
if [ -z "$1" ]; then
    # No arguments, check stdin
    if [ ! -t 0 ]; then
        SETUP_FILE="setup.yml"
        INPUT=$(cat)
    else
        echo "Error: No input provided"
        echo "Usage: $0 <PR number>            e.g. $0 984"
        echo "   or: $0 ci_output.txt [setup.yml]"
        echo "   or: $0 setup.yml ci_output.txt"
        echo "   or: $0 < ci_output.txt"
        echo "   or: $0 setup.yml \"CI output text\""
        exit 1
    fi
elif [ -f "$1" ] && grep -q "Variable:.*current_value:" "$1" 2>/dev/null; then
    # First arg is a CI output file
    INPUT=$(cat "$1")
    SETUP_FILE="${2:-setup.yml}"
elif [ -f "$1" ]; then
    # First arg is setup file
    SETUP_FILE="$1"
    if [ -n "$2" ]; then
        if [ -f "$2" ]; then
            INPUT=$(cat "$2")
        else
            INPUT="$2"
        fi
    elif [ ! -t 0 ]; then
        INPUT=$(cat)
    else
        echo "Error: No CI output provided"
        echo "Usage: $0 setup.yml ci_output.txt"
        echo "   or: $0 setup.yml \"CI output text\""
        echo "   or: cat ci_output.txt | $0 setup.yml"
        exit 1
    fi
else
    echo "Error: File $1 not found"
    exit 1
fi

# Check if setup.yml exists
if [ ! -f "$SETUP_FILE" ]; then
    echo "Error: $SETUP_FILE not found"
    echo "Available setup directories:"
    find . -maxdepth 2 -name "setup.yml" -type f 2>/dev/null | sed 's|^\./||'
    exit 1
fi

# Parse and update values
updates_keys=()
updates_vals=()
count=0

while IFS= read -r line; do
    # Match lines like: "Variable: a_ice, current_value: 0.19109711397780874, master_value: ..."
    if [[ "$line" =~ Variable:\ ([^,]+),\ current_value:\ ([^,]+), ]]; then
        var_name="${BASH_REMATCH[1]}"
        var_name=$(echo "$var_name" | xargs)  # trim whitespace
        current_val="${BASH_REMATCH[2]}"
        current_val=$(echo "$current_val" | xargs)  # trim whitespace

        updates_keys+=("$var_name")
        updates_vals+=("$current_val")
        ((count++))
    fi
done <<< "$INPUT"

if [ $count -eq 0 ]; then
    echo "Error: No variable updates found in input"
    echo "Expected format: Variable: <name>, current_value: <value>, master_value: ..."
    exit 1
fi

echo "Found $count variable(s) to update:"
for i in "${!updates_keys[@]}"; do
    echo "  ${updates_keys[$i]}: ${updates_vals[$i]}"
done
echo ""

# Apply updates to setup.yml
uname_s=$(uname -s)
for i in "${!updates_keys[@]}"; do
    var="${updates_keys[$i]}"
    value="${updates_vals[$i]}"
    # Escape special characters in value for sed
    escaped_value=$(echo "$value" | sed 's/[\/&]/\\&/g')

    if [ "$uname_s" = "Darwin" ]; then
        sed -i '' -e "s/^\\([ ]*${var}:\\).*/\\1 ${escaped_value}/" "$SETUP_FILE"
    else
        sed -i -e "s/^\\([ ]*${var}:\\).*/\\1 ${escaped_value}/" "$SETUP_FILE"
    fi
done

echo "✓ Truth values updated successfully in $SETUP_FILE"
