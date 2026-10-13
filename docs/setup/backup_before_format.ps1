<#
.SYNOPSIS
  Backup of everything needed to rebuild this PC (Windows 11 reinstall) and restore the
  pbl-analysis setup exactly. Companion of docs/setup/pc_rebuild_checklist.md.

.USAGE
  Run in a normal (non-admin) PowerShell, with an external drive connected, e.g. E:\
    powershell -ExecutionPolicy Bypass -File C:\WORK\dev\pbl-analysis\docs\setup\backup_before_format.ps1 -BackupRoot E:\
  Optional switches:
    -IncludeUserFolders     also copy Desktop, Documents, Downloads
    -IncludeClaudeAppData   also copy the whole Claude Desktop app-data folder (can be large)

.NOTES
  Read-only on the source: it only copies. Writes everything to <BackupRoot>\PC_BACKUP_<timestamp>\
  Check backup_log.txt and REPO_STATUS.txt at the end.
#>
param(
  [Parameter(Mandatory = $true)][string]$BackupRoot,
  [switch]$IncludeUserFolders,
  [switch]$IncludeClaudeAppData
)

$ErrorActionPreference = "Continue"
$stamp = Get-Date -Format "yyyyMMdd_HHmm"
$dest  = Join-Path $BackupRoot "PC_BACKUP_$stamp"
$cfg   = Join-Path $dest "config"
New-Item -ItemType Directory -Force -Path $dest, $cfg | Out-Null
$log = Join-Path $dest "backup_log.txt"

function Log($msg) { $line = "[{0}] {1}" -f (Get-Date -Format "HH:mm:ss"), $msg; Write-Host $line; Add-Content -Path $log -Value $line }
function Has($cmd) { return [bool](Get-Command $cmd -ErrorAction SilentlyContinue) }
function Copy-Tree($src, $dst, $extra = @()) {
  if (-not (Test-Path $src)) { Log "  skip (not found): $src"; return }
  $rcArgs = @($src, $dst, "/E", "/R:1", "/W:1", "/NFL", "/NDL", "/NP", "/LOG+:$log") + $extra
  & robocopy @rcArgs | Out-Null
  if ($LASTEXITCODE -ge 8) { Log "  ERROR robocopy ($LASTEXITCODE): $src" } else { Log "  copied: $src" }
}

Log "Backup destination: $dest"

# 1. Git repositories: uncommitted changes and unpushed commits -------------------------------
Log "1. Checking git repositories under C:\WORK\dev"
$repoReport = Join-Path $dest "REPO_STATUS.txt"
"Repository status at $stamp (empty sections = clean)" | Set-Content $repoReport
if (Has git) {
  Get-ChildItem C:\WORK\dev -Directory -ErrorAction SilentlyContinue | Where-Object { Test-Path (Join-Path $_.FullName ".git") } | ForEach-Object {
    $r = $_.FullName
    Add-Content $repoReport "`n=== $r"
    Add-Content $repoReport "--- remote:";             git -C $r remote -v | Add-Content $repoReport
    Add-Content $repoReport "--- uncommitted changes:"; git -C $r status --short | Add-Content $repoReport
    Add-Content $repoReport "--- commits not on any remote:"; git -C $r log --branches --not --remotes --oneline | Add-Content $repoReport
    Add-Content $repoReport "--- branches:";           git -C $r branch -a | Add-Content $repoReport
  }
  Log "  written: REPO_STATUS.txt (CHECK IT: anything listed must be committed/pushed or is only in this backup)"
}

# 2. Software inventory -----------------------------------------------------------------------
Log "2. Software inventory"
if (Has winget) {
  winget export -o (Join-Path $cfg "winget_packages.json") --accept-source-agreements | Out-Null
  winget list --accept-source-agreements | Out-File (Join-Path $cfg "winget_list.txt") -Encoding utf8
  Log "  winget export + list done"
}
$keys = "HKLM:\Software\Microsoft\Windows\CurrentVersion\Uninstall\*",
        "HKLM:\Software\WOW6432Node\Microsoft\Windows\CurrentVersion\Uninstall\*",
        "HKCU:\Software\Microsoft\Windows\CurrentVersion\Uninstall\*"
Get-ItemProperty $keys -ErrorAction SilentlyContinue | Where-Object DisplayName |
  Select-Object DisplayName, DisplayVersion, Publisher, InstallLocation |
  Sort-Object DisplayName | Export-Csv (Join-Path $cfg "installed_programs.csv") -NoTypeInformation -Encoding utf8
Log "  installed_programs.csv done"

# 3. Tool versions and install locations ------------------------------------------------------
Log "3. Tool versions"
$v = Join-Path $cfg "versions.txt"
"Tool versions at $stamp" | Set-Content $v
foreach ($c in "git", "gh", "claude", "code", "Rscript", "R") {
  $cmd = Get-Command $c -ErrorAction SilentlyContinue
  if ($cmd) { Add-Content $v "`n$c -> $($cmd.Source)" } else { Add-Content $v "`n$c -> NOT FOUND" }
}
if (Has git)     { git --version 2>&1 | Add-Content $v }
if (Has gh)      { gh --version 2>&1 | Add-Content $v }
if (Has claude)  { claude --version 2>&1 | Add-Content $v }
if (Has code)    { code --version 2>&1 | Add-Content $v }
if (Has Rscript) { Rscript --version 2>&1 | Add-Content $v }
Add-Content $v "`nR installations:"; Get-ChildItem "C:\Program Files\R" -Directory -ErrorAction SilentlyContinue | ForEach-Object { Add-Content $v "  $($_.FullName)" }
Add-Content $v "Rtools installations:"; Get-ChildItem "C:\" -Directory -Filter "rtools*" -ErrorAction SilentlyContinue | ForEach-Object { Add-Content $v "  $($_.FullName)" }
Add-Content $v "`nWindows:"; (Get-CimInstance Win32_OperatingSystem | Select-Object Caption, Version | Out-String) | Add-Content $v
Log "  versions.txt done"

# 4. Environment variables (PATH) -------------------------------------------------------------
Log "4. Environment variables"
[Environment]::GetEnvironmentVariable("Path", "User")    | Set-Content (Join-Path $cfg "PATH_user.txt")
[Environment]::GetEnvironmentVariable("Path", "Machine") | Set-Content (Join-Path $cfg "PATH_machine.txt")
Get-ChildItem Env: | Where-Object { $_.Name -match '^R_|RENV|CONFIG|GITHUB|CLAUDE' -and $_.Name -notmatch 'TOKEN|KEY|SECRET' } |
  Format-Table -AutoSize | Out-String | Set-Content (Join-Path $cfg "env_vars_selected.txt")

# 5. Git and GitHub CLI configuration -------------------------------------------------------------
Log "5. Git config"
if (Test-Path "$HOME\.gitconfig") { Copy-Item "$HOME\.gitconfig" (Join-Path $cfg "gitconfig") }
if (Has git) { git config --global --list | Set-Content (Join-Path $cfg "git_config_global.txt") }
if (Has gh)  { gh auth status 2>&1 | Set-Content (Join-Path $cfg "gh_auth_status.txt") }   # token is masked by gh

# 6. VS Code -------------------------------------------------------------------------------------
Log "6. VS Code settings and extensions"
$vs = Join-Path $cfg "vscode"; New-Item -ItemType Directory -Force -Path $vs | Out-Null
if (Has code) { code --list-extensions --show-versions | Set-Content (Join-Path $vs "extensions.txt") }
foreach ($f in "settings.json", "keybindings.json") {
  $p = Join-Path $env:APPDATA "Code\User\$f"; if (Test-Path $p) { Copy-Item $p $vs }
}
Copy-Tree (Join-Path $env:APPDATA "Code\User\snippets") (Join-Path $vs "snippets")

# 7. Claude Code (CLI / VS Code extension) settings ----------------------------------------------
Log "7. Claude Code settings"
Copy-Tree "$HOME\.claude" (Join-Path $cfg "dot_claude")
if (Test-Path "$HOME\.claude.json") { Copy-Item "$HOME\.claude.json" (Join-Path $cfg "dot_claude.json") }

# 8. Claude Desktop ------------------------------------------------------------------------------
Log "8. Claude Desktop config"
$cd = Join-Path $cfg "claude_desktop"; New-Item -ItemType Directory -Force -Path $cd | Out-Null
Get-ChildItem (Join-Path $env:APPDATA "Claude") -File -Filter *.json -ErrorAction SilentlyContinue | Copy-Item -Destination $cd
if ($IncludeClaudeAppData) {
  Copy-Tree (Join-Path $env:APPDATA "Claude") (Join-Path $dest "AppData_Claude") @("/XD", "Cache", "Code Cache", "GPUCache", "DawnCache", "blob_storage")
}

# 9. R: globally installed packages (all R versions found) ----------------------------------------
Log "9. R package lists"
Get-ChildItem "C:\Program Files\R" -Directory -ErrorAction SilentlyContinue | ForEach-Object {
  $rs = Join-Path $_.FullName "bin\Rscript.exe"
  if (Test-Path $rs) {
    $out = (Join-Path $cfg ("R_packages_" + $_.Name + ".csv")) -replace '\\', '/'
    & $rs -e "ip <- installed.packages()[, c('Package','Version','LibPath')]; write.csv(ip, '$out', row.names = FALSE)" 2>&1 | Out-Null
    Log "  $($_.Name) package list"
  }
}
foreach ($f in "$HOME\.Rprofile", "$HOME\.Renviron", "$HOME\Documents\.Rprofile", "$HOME\Documents\.Renviron") {
  if (Test-Path $f) { Copy-Item $f (Join-Path $cfg ("R_" + (Split-Path $f -Leaf))) ; Log "  copied $f (may contain secrets: keep the backup private)" }
}

# 10. Remote-access tools (cluster) ---------------------------------------------------------------
Log "10. SSH / WinSCP / MobaXterm / PuTTY settings"
Copy-Tree "$HOME\.ssh" (Join-Path $cfg "dot_ssh")
if (Test-Path "HKCU:\Software\Martin Prikryl") { reg export "HKCU\Software\Martin Prikryl" (Join-Path $cfg "winscp_sessions.reg") /y | Out-Null; Log "  WinSCP sessions exported" }
if (Test-Path "HKCU:\Software\SimonTatham")    { reg export "HKCU\Software\SimonTatham" (Join-Path $cfg "putty_sessions.reg") /y | Out-Null; Log "  PuTTY sessions exported" }
Copy-Tree (Join-Path $env:APPDATA "MobaXterm") (Join-Path $cfg "MobaXterm_AppData")
Copy-Tree "$HOME\Documents\MobaXterm" (Join-Path $cfg "MobaXterm_Documents")

# 11. C:\WORK (all projects, data, secrets) -------------------------------------------------------
Log "11. Copying C:\WORK (excluding rebuildable renv libraries)"
$xd = @(Get-ChildItem C:\WORK -Recurse -Directory -ErrorAction SilentlyContinue |
        Where-Object { $_.Parent -and $_.Parent.Name -eq "renv" -and $_.Name -in @("library", "staging", "sandbox") } |
        ForEach-Object { $_.FullName })
Log "  excluded renv folders: $($xd.Count)"
$extra = @(); if ($xd.Count -gt 0) { $extra = @("/XD") + $xd }
Copy-Tree "C:\WORK" (Join-Path $dest "WORK") $extra

# 12. Optional user folders ------------------------------------------------------------------------
if ($IncludeUserFolders) {
  Log "12. User folders"
  foreach ($f in "Desktop", "Documents", "Downloads") { Copy-Tree "$HOME\$f" (Join-Path $dest "User_$f") }
}

# Summary --------------------------------------------------------------------------------------------
$size = (Get-ChildItem $dest -Recurse -File -ErrorAction SilentlyContinue | Measure-Object Length -Sum).Sum / 1GB
Log ("DONE. Backup size: {0:N2} GB in {1}" -f $size, $dest)
Log "Next: open REPO_STATUS.txt and backup_log.txt; verify that WORK\dev\pbl-analysis\data\raw and WORK\secrets are in the backup."
