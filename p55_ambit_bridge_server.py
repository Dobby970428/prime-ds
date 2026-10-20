"""
PRIME-DS local AMBIT2 REST bridge.

Implements the PRIME-DS custom endpoint expected by ambit_client.py:
  POST /ambit2/query/property

Uses the same RDKit (Q)SAR proxy as offline fallback so results stay consistent.
Optional — PRIME-DS works without this server (direct RDKit proxy in ambit_client).
"""

# 로컬 흉내내기 서버(별도 실행)
# run_ambit_bridge.bat을 켜면 화면에 'Connected(AMBIT2)'로 뜨지만, **그 서버는 진짜 OECD Toolbox가 아니라 p55가 만든 흉내내기 서버**임.
# 내부적으로 RDKit 계산 결과를 AMBIT 형식으로 포장해 돌려줌.
# 어느 경로로 가든 실제 계산은 RDKit임.

from __future__ import annotations

import json
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any

from p53_oecd_rdkit_fallback import rdkit_oecd_qsar_proxy
from p01_prime_config import get_ambit_bridge_bind, get_ambit_bridge_port

_BRIDGE_VERSION = "1.0.0"


def proxy_to_ambit_api_response(proxy: dict[str, str]) -> dict[str, Any]:
    """Map internal RDKit proxy keys → AMBIT2 REST JSON (ambit_client contract)."""
    return {
        "hydrolysis_rate": proxy.get("Hydrolysis_Stability", "N/A"),
        "biodegradability": proxy.get("Biodegradability", "N/A"),
        "mutagenicity_alert": proxy.get("OECD_Alert_Found", "None detected"),
        "source": "PRIME-DS-AMBIT-Bridge",
        "engine": "PRIME-DS RDKit (Q)SAR proxy (local AMBIT2-compatible bridge)",
        "bridge_version": _BRIDGE_VERSION,
    }


class AmbitBridgeHandler(BaseHTTPRequestHandler):
    server_version = f"PRIME-DS-AMBIT-Bridge/{_BRIDGE_VERSION}"

    def log_message(self, fmt: str, *args) -> None:
        sys.stderr.write("%s - [%s] %s\n" % (self.address_string(), self.log_date_time_string(), fmt % args))

    def _send_json(self, status: int, payload: dict[str, Any]) -> None:
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:
        path = self.path.split("?", 1)[0]
        if path in ("/ambit2", "/ambit2/", "/ambit2/health"):
            self._send_json(
                200,
                {
                    "status": "ok",
                    "service": "PRIME-DS AMBIT2 bridge",
                    "version": _BRIDGE_VERSION,
                    "endpoint": "/ambit2/query/property",
                },
            )
            return
        self._send_json(404, {"error": "Not found"})

    def do_POST(self) -> None:
        path = self.path.split("?", 1)[0]
        if path != "/ambit2/query/property":
            self._send_json(404, {"error": "Not found"})
            return

        length = int(self.headers.get("Content-Length", 0) or 0)
        raw = self.rfile.read(length).decode("utf-8", errors="replace") if length else "{}"
        try:
            payload = json.loads(raw)
        except json.JSONDecodeError:
            self._send_json(400, {"error": "Invalid JSON body"})
            return

        smiles = str(payload.get("smiles") or "").strip()
        if not smiles:
            self._send_json(400, {"error": "Missing 'smiles' in request body"})
            return

        proxy = rdkit_oecd_qsar_proxy(smiles)
        self._send_json(200, proxy_to_ambit_api_response(proxy))


def run_server(*, host: str | None = None, port: int | None = None) -> None:
    bind = host or get_ambit_bridge_bind()
    listen_port = int(port if port is not None else get_ambit_bridge_port())
    httpd = ThreadingHTTPServer((bind, listen_port), AmbitBridgeHandler)
    base = f"http://{bind}:{listen_port}/ambit2"
    print("=" * 60)
    print("  PRIME-DS AMBIT2 Bridge Server")
    print(f"  Health : {base}/health")
    print(f"  Query  : POST {base}/query/property")
    print("  Stop   : Ctrl+C")
    print("=" * 60)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nBridge server stopped.")
        httpd.server_close()


if __name__ == "__main__":
    run_server()
