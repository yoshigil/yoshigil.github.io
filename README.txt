UPDATING johangilchrist.com — GIT WORKFLOW
===========================================
This site is a static HTML site hosted on GitHub Pages.
Repo:   https://github.com/yoshigil/yoshigil.github.io
Branch: main
Domain: johangilchrist.com (set via the CNAME file — don't delete it)

There is no build step. Whatever HTML/CSS/JS/image files are on the
`main` branch are what's live. GitHub Pages typically republishes the
site within 1-2 minutes of a push (occasionally up to ~10 min).


BASIC WORKFLOW (every time you make changes)
----------------------------------------------
1. Open Terminal and go to the site folder:
       cd path/to/Personal_Website

2. See what you've changed:
       git status
   (Modified files show as "modified", new files as "untracked",
   deleted files as "deleted".)

3. (Optional) Review the actual changes before committing:
       git diff

4. Stage the files you want to include in this update:
       git add <file1> <file2> ...
   or, to stage everything you've changed:
       git add -A

5. Commit the staged changes with a short description:
       git commit -m "Describe what changed, e.g. Add vortex ring model page"

6. Push to GitHub (this is what actually publishes the update):
       git push

That's it — GitHub Pages rebuilds automatically after step 6. No
separate deploy command needed.


CHECKING THE LIVE SITE
-----------------------
After pushing, wait a minute or two, then reload
https://johangilchrist.com (hard-refresh with Cmd+Shift+R if you don't
see the change, since browsers cache aggressively).

You can also confirm the push succeeded and see build status at:
https://github.com/yoshigil/yoshigil.github.io/actions


UNDOING MISTAKES
------------------
- Unstage a file you `git add`ed by accident:
      git restore --staged <file>

- Discard local edits to a file (revert to last commit):
      git restore <file>

- Fix the message or add a forgotten file to your LAST commit
  (only if you haven't pushed yet):
      git add <forgotten-file>
      git commit --amend

- Already pushed a bad commit? Don't rewrite history on a shared repo
  — instead make a new commit that fixes it:
      git commit -m "Fix typo from previous commit"
      git push


TYPICAL EXAMPLE
------------------
    cd path/to/Personal_Website
    git status
    git add research/Volcanology/volcanology.html research.html
    git commit -m "Fix broken links after moving volcanology page"
    git push


TROUBLESHOOTING
------------------
"fatal: Unable to create '.git/index.lock': File exists."
    An earlier git process was interrupted and left a stale lock file.
    Quit any open Terminal/GUI git clients, then delete it manually:
        rm .git/index.lock
    (If `rm` also refuses with a permissions error, delete the file
    /path/to/Personal_Website/.git/index.lock in Finder, or restart
    your Mac if it's locked by a stuck process.)

"Your branch is behind 'origin/main'"
    Someone (or another device) pushed changes you don't have locally.
    Pull them first, then push:
        git pull
        git push

Nothing shows up on the live site after pushing
    - Confirm the push actually went through: `git log --oneline -1`
      should match the latest commit shown on GitHub.
    - Hard-refresh the browser (Cmd+Shift+R) — GitHub Pages / your
      browser cache is the most common culprit.
    - Check https://github.com/yoshigil/yoshigil.github.io/actions
      for a failed Pages build.
