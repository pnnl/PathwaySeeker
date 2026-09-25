"""Install the PathwaySeeker skill for Claude Code and Codex, and register the MCP server with Codex.

    pathwayseeker setup            # every assistant found on this machine
    pathwayseeker setup claude     # just Claude Code
    pathwayseeker setup codex      # just Codex
"""

import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional

SKILL_SOURCE = Path(__file__).parent / "skill" / "SKILL.md"
CODEX_MCP_BLOCK = """
[mcp_servers.pathwayseeker]
command = "{exe}"
args = ["mcp"]
"""


def pathwayseeker_exe() -> str:
    """Absolute path of the pathwayseeker command in this environment."""
    found = shutil.which("pathwayseeker")
    if found:
        return str(Path(found).resolve())
    return str(Path(sys.executable).parent / "pathwayseeker")


def _skill_text(exe: str) -> str:
    text = SKILL_SOURCE.read_text()
    note = (f"\n> On this machine the `pathwayseeker` command is `{exe}`. If `pathwayseeker` is not "
            f"found in the shell, use that full path.\n")
    head, sep, body = text.partition("\n# PathwaySeeker\n")
    return head + sep + note + body if sep else text + note


def install_skill(skills_dir: Path, exe: str) -> Path:
    dest = skills_dir / "pathwayseeker" / "SKILL.md"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(_skill_text(exe))
    return dest


def register_codex_mcp(codex_home: Path, exe: str) -> str:
    config = codex_home / "config.toml"
    existing = config.read_text() if config.exists() else ""
    if "[mcp_servers.pathwayseeker]" in existing:
        return f"already registered in {config}"
    config.parent.mkdir(parents=True, exist_ok=True)
    with open(config, "a") as f:
        f.write(CODEX_MCP_BLOCK.format(exe=exe))
    return f"added to {config}"


def setup(targets: Optional[List[str]] = None, home: Optional[Path] = None) -> dict:
    home = Path(home) if home else Path.home()
    exe = pathwayseeker_exe()
    claude_dir, codex_dir = home / ".claude", home / ".codex"
    if not targets:
        targets = [t for t, d in (("claude", claude_dir), ("codex", codex_dir)) if d.exists()]
        if not targets:
            targets = ["claude", "codex"]
    done = {"command": exe}
    if "claude" in targets:
        done["claude"] = {"skill": str(install_skill(claude_dir / "skills", exe)),
                          "next": "Start a new Claude Code session and ask a question about your data."}
    if "codex" in targets:
        done["codex"] = {"skill": str(install_skill(codex_dir / "skills", exe)),
                         "mcp_server": register_codex_mcp(codex_dir, exe),
                         "next": "Start a new Codex session and ask a question about your data."}
    try:
        subprocess.run([exe, "--help"], capture_output=True, check=True, timeout=60)
    except Exception as e:  # the assistant needs to be able to run the command
        done["warning"] = f"Could not run {exe}: {e}"
    return done
