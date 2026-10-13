# PC rebuild checklist (format + Windows 11)

Goal: after formatting, get the PC back to the exact working state of 2026-10-13 — software, configuration,
and the `pbl-analysis` project — and verify it.
This file lives in the GitHub repo, so it survives the format: open it from
https://github.com/Giacomo-Nicolini/pbl-analysis/blob/main/docs/setup/pc_rebuild_checklist.md

Tick the boxes as you go. Commands are for PowerShell unless stated otherwise.

---

## Part A — Before formatting (on the current PC)

### A1. Make sure the project is fully on GitHub
- [ ] In `C:\WORK\dev\pbl-analysis`: `git status` shows nothing to commit; `git push` says "Everything up-to-date".
- [ ] Optional: push the old spec branches too (they are already merged into `main`, so this is only for history):
      `git push origin spec-01-readers spec-02-pblh-model`

### A2. Run the backup script (external drive, e.g. `E:\`)
- [ ] Connect an external drive with enough space (C:\WORK contains the old R projects and the raw data).
- [ ] Run (normal PowerShell, not admin):
      ```powershell
      powershell -ExecutionPolicy Bypass -File C:\WORK\dev\pbl-analysis\docs\setup\backup_before_format.ps1 -BackupRoot E:\ -IncludeUserFolders
      ```
      Add `-IncludeClaudeAppData` to also keep the Claude Desktop app data (local Cowork sessions); it can be large.
- [ ] Open `E:\PC_BACKUP_<date>\REPO_STATUS.txt`: every repo must show **no** uncommitted changes and **no** unpushed
      commits (otherwise commit/push and rerun).
- [ ] Open `backup_log.txt`: no `ERROR` lines.

What the script saves (in `E:\PC_BACKUP_<date>\`):

| Item | Where in the backup | Why it matters |
|---|---|---|
| All of `C:\WORK` (old R projects, `dev`, `secrets`) | `WORK\` | **Not in git**: `data/raw` (ICOS zips, ERA5, target-area KML), `C:\WORK\R`, `C:\WORK\secrets` (ERA5 API key), EBC ancillary files |
| Installed software list | `config\winget_packages.json`, `installed_programs.csv` | reinstall the same programs |
| Tool versions and paths | `config\versions.txt` | e.g. how Claude Code was installed (path of `claude`) |
| User/Machine PATH | `config\PATH_*.txt` | R 4.4.1 on PATH |
| Git config, gh status | `config\gitconfig`, `git_config_global.txt` | identity, `init.defaultBranch`, `core.autocrlf` |
| VS Code extensions, settings | `config\vscode\` | same editor setup |
| Claude Code settings | `config\dot_claude\`, `dot_claude.json` | permissions, preferences |
| Claude Desktop config | `config\claude_desktop\` | app settings |
| R package lists per R version | `config\R_packages_R-*.csv` | old projects outside renv |
| SSH / WinSCP / PuTTY / MobaXterm sessions | `config\dot_ssh`, `*.reg`, `MobaXterm_*` | cluster access (`dibaf-fe2`) |

### A3. Things the script cannot save — check manually
- [ ] Browser bookmarks/passwords (sign in to browser sync, or export).
- [ ] Mendeley: library is in your account (sync); local PDF folder if you keep one outside C:\WORK.
- [ ] VPN / CMCC / university clients and their settings; licence keys (Office, etc.).
- [ ] Anything stored outside `C:\WORK` and your user folders (e.g. other drives' folders you use).
- [ ] Write down the cluster host name/IP you use in WinSCP/MobaXterm (also in the exported sessions).
- [ ] Windows 11 compatibility: run Microsoft "PC Health Check" (TPM 2.0, Secure Boot, CPU).
- [ ] Copy the backup a second time (cloud or second drive) if possible. Only then format.

### A4. Reference values to verify the restore (current state)
- `testthat::test_dir("tests/testthat")`: **375 tests, 0 failures**.
- `scripts/02_run_pblh_model.R` on BE-Bra, default case: convective median **704 m** (p5–p95 72–1609),
  stable median **271 m**, all rows median **423 m**; 100,601 rows valid in all cases.

---

## Part B — After installing Windows 11

### B1. Base system
- [ ] Windows Update until nothing is left; reboot.
- [ ] Check `winget --version` (included in Windows 11). If `msstore` certificate errors appear, always add
      `--source winget` to the commands below (as on Windows 10).

### B2. Core tools (same as before)
```powershell
winget install --id Git.Git -e --source winget
winget install --id GitHub.cli -e --source winget
winget install --id Microsoft.VisualStudioCode -e --source winget
winget install --id Anthropic.Claude -e --source winget        # Claude Desktop; if the id is not found: https://claude.ai/download
```
- [ ] Close and reopen PowerShell. Check: `git --version`, `gh --version`, `code --version`.
- [ ] Other programs: compare with `config\winget_packages.json` / `installed_programs.csv`.
      You can reinstall many at once with `winget import -i E:\PC_BACKUP_<date>\config\winget_packages.json`
      (review the list first; it installs everything it recognises).

### B3. Claude Code (CLI, used by VS Code)
- [ ] Look at `config\versions.txt` → line `claude -> <path>` to see how it was installed:
  - path under `...\.local\bin\` (native installer) → in PowerShell: `irm https://claude.ai/install.ps1 | iex`
  - path under `...\npm\` (npm) → install Node.js LTS, then `npm install -g @anthropic-ai/claude-code`
  - If in doubt, follow the current instructions at https://docs.claude.com (Claude Code → Setup).
- [ ] `claude --version` works; run `claude` once and sign in with your Claude account.
- [ ] Optional: copy back `config\dot_claude\settings.json` to `%USERPROFILE%\.claude\` (only settings files;
      do not copy old session folders).

### B4. Git and GitHub
```powershell
git config --global user.name "Giacomo Nicolini"
git config --global user.email "giacomo.nicolini@cmcc.it"
git config --global init.defaultBranch main
git config --global core.autocrlf true
gh auth login        # GitHub.com → HTTPS → Yes → Login with a web browser
gh auth status
```

### B5. R toolchain (must match the cluster: R 4.4.1)
- [ ] R 4.4.1: https://cran.r-project.org/bin/windows/base/old/4.4.1/R-4.4.1-win.exe (default folder).
- [ ] Rtools44: https://cran.r-project.org/bin/windows/Rtools/rtools44/rtools.html (default `C:\rtools44`).
- [ ] Optional, only for old projects: the other R versions listed in `config\versions.txt` (4.6.1 + Rtools45, …).
- [ ] Put R 4.4.1 on the user PATH:
      ```powershell
      $rbin = "C:\Program Files\R\R-4.4.1\bin"
      $old  = [Environment]::GetEnvironmentVariable("Path", "User")
      if ($old -notlike "*$rbin*") { [Environment]::SetEnvironmentVariable("Path", "$old;$rbin", "User") }
      ```
- [ ] Reopen PowerShell: `Rscript -e "R.version.string; Sys.which('make')"` → R 4.4.1 and `C:\rtools44\usr\bin\make.exe`.
- [ ] Optional: RStudio, then Tools → Global Options → General → R version → 4.4.1.

### B6. Restore C:\WORK
- [ ] Copy `E:\PC_BACKUP_<date>\WORK` to `C:\WORK` (this restores `C:\WORK\R`, `C:\WORK\secrets`, and
      `C:\WORK\dev\...` including `pbl-analysis\data\raw`).
      ```powershell
      robocopy E:\PC_BACKUP_<date>\WORK C:\WORK /E
      ```
- [ ] Alternative for the project (clean clone) — then copy only the gitignored data back:
      ```powershell
      gh repo clone Giacomo-Nicolini/pbl-analysis C:\WORK\dev\pbl-analysis
      robocopy E:\PC_BACKUP_<date>\WORK\dev\pbl-analysis\data C:\WORK\dev\pbl-analysis\data /E
      ```
- [ ] In `C:\WORK\dev\pbl-analysis`: `git status` clean, `git pull` up to date.

### B7. Rebuild the project R environment (renv)
```powershell
cd C:\WORK\dev\pbl-analysis
Rscript -e "renv::restore(prompt = FALSE)"
Rscript -e "renv::install(c('languageserver','httpgd'))"     # editor helpers, not recorded in renv.lock
Rscript -e "renv::status()"                                 # expected: consistent state
```
(renv bootstraps itself from `renv/activate.R`; packages come from Posit Package Manager as set in `.Rprofile`.)

### B8. VS Code
- [ ] Extensions: **R** (REditorSupport, id `REditorSupport.r`) and **Claude Code** (Anthropic, id `anthropic.claude-code`),
      plus any other in `config\vscode\extensions.txt`:
      ```powershell
      Get-Content E:\PC_BACKUP_<date>\config\vscode\extensions.txt | ForEach-Object { code --install-extension ($_ -split '@')[0] }
      ```
- [ ] Copy back `settings.json` / `keybindings.json` from `config\vscode\` to `%APPDATA%\Code\User\` (or use VS Code
      Settings Sync with your GitHub account).
- [ ] `cd C:\WORK\dev\pbl-analysis; code .` → Claude Code panel opens and reads `CLAUDE.md`.

### B9. Cluster access
- [ ] Reinstall WinSCP / MobaXterm / PuTTY as before; import sessions (`*.reg` double-click, or MobaXterm ini).
- [ ] Copy `config\dot_ssh` back to `%USERPROFILE%\.ssh` if you used SSH keys.
- [ ] Test: connect to the cluster and list `/bigdata/g.nicolini/ETC/L2_DATA_backup/L2I_2026/`.

### B10. Claude Desktop (Cowork)
- [ ] Sign in with the same account. Choose the model in the model selector.
- [ ] For the PBLH work, start a Cowork session and grant access to `C:\WORK\dev\pbl-analysis`.
      The project memory lives in the repo: ask Claude to read `docs/STATUS.md`, `docs/project_plan.md` and `CLAUDE.md`
      first — everything decided so far is there.
- [ ] If you backed up the app data (`-IncludeClaudeAppData`), keep it as an archive; restoring old local sessions is
      not needed to continue the project.

### B11. Verify (compare with A4)
```powershell
cd C:\WORK\dev\pbl-analysis
Rscript -e "testthat::test_dir('tests/testthat')"     # expect 375 tests, 0 failures
Rscript scripts/01_check_inputs.R                        # BE-Bra input report, no errors
Rscript scripts/02_run_pblh_model.R                      # default convective median 704 m, stable 271 m
```
- [ ] Numbers match A4 → restore complete. If not, check R version (`4.4.1`), `renv::status()`, and that
      `data\raw` was restored completely.

### B12. Clean-up
- [ ] Keep the external backup until the first successful cluster run after the rebuild.
- [ ] Delete `C:\WORK\secrets` copies from any place that is not private.
