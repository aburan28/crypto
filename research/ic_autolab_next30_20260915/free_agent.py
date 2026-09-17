"""Harbor adapter for the preinstalled, credential-free OpenCode loop."""
import shlex
from harbor.agents.installed.opencode import OpenCode

class FreeOpenCode(OpenCode):
    @staticmethod
    def name():
        return "ic-free-opencode"

    def get_version_command(self):
        return "/usr/local/bin/opencode --version"

    async def install(self, environment):
        await self.exec_as_root(environment, command="mkdir -p /logs/agent; chown researcher:researcher /logs/agent; chmod 755 /logs/verifier")
        await self.exec_as_agent(environment, command="/usr/local/bin/opencode --version")

    async def run(self, instruction, environment, context):
        await self.exec_as_agent(environment,
            command="python3 /opt/harness/agent_loop.py --seconds 3600 --instruction " + shlex.quote(instruction),
            cwd="/app")
