import http.client
import json
import threading

import pytest

from minimol.server import make_server


@pytest.fixture(scope="module")
def server(tmp_path_factory):
    path = tmp_path_factory.mktemp("files") / "water.xyz"
    path.write_text("3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n")
    srv = make_server("127.0.0.1", 0, str(path))
    thread = threading.Thread(target=srv.serve_forever, daemon=True)
    thread.start()
    yield srv
    srv.shutdown()
    srv.server_close()


def request(server, method, path, body=None, headers=None, host=None):
    port = server.server_address[1]
    conn = http.client.HTTPConnection("127.0.0.1", port, timeout=60)
    h = {"Host": host or f"127.0.0.1:{port}"}
    if body is not None:
        h["Content-Type"] = "application/json"
        h["X-Minimol"] = "1"
    h.update(headers or {})
    conn.request(method, path, body=json.dumps(body) if body is not None else None, headers=h)
    res = conn.getresponse()
    data = res.read()
    conn.close()
    return res, data


def api(server, path, body):
    res, data = request(server, "POST", path, body)
    return res.status, json.loads(data)


def test_index_and_static_files(server):
    res, data = request(server, "GET", "/")
    assert res.status == 200 and b"minimol" in data
    assert "script-src 'self'" in res.getheader("Content-Security-Policy")
    res, _ = request(server, "GET", "/static/js/app.js")
    assert res.status == 200 and res.getheader("Content-Type").startswith("text/javascript")
    res, _ = request(server, "GET", "/static/vendor/three.module.js")
    assert res.status == 200


@pytest.mark.parametrize("path", ["/static/../server.py", "/static/%2e%2e/server.py", "/etc/passwd", "/static/js/../../server.py"])
def test_no_path_traversal(server, path):
    res, _ = request(server, "GET", path)
    assert res.status == 404


def test_foreign_host_is_rejected(server):
    res, _ = request(server, "GET", "/", host="evil.example:80")
    assert res.status == 403
    res, _ = request(server, "POST", "/api/info", {}, host="evil.example")
    assert res.status == 403


def test_api_requires_custom_header(server):
    res, _ = request(server, "POST", "/api/info", {}, headers={"X-Minimol": "0"})
    assert res.status == 403


def test_info_lists_methods(server):
    status, data = api(server, "/api/info", {})
    assert status == 200 and data["ok"]
    assert data["methods"][:3] == ["MMFF94", "MMFF94s", "UFF"]


def test_initial_file(server):
    status, data = api(server, "/api/initial", {})
    assert status == 200 and data["name"] == "water.xyz"
    assert [a["el"] for a in data["molecule"]["atoms"]] == ["O", "H", "H"]


def test_smiles_addh_optimize_export_flow(server):
    status, data = api(server, "/api/smiles", {"smiles": "CCO"})
    assert status == 200
    mol = data["molecule"]
    status, opt = api(server, "/api/optimize", {"molecule": mol, "method": "MMFF94"})
    assert status == 200 and opt["converged"] and opt["method"] == "MMFF94"
    status, out = api(server, "/api/export", {"molecule": opt["molecule"], "format": "xyz", "title": "eth\nanol"})
    assert status == 200 and out["text"].splitlines()[:2] == ["9", "eth anol"]


def test_user_errors_are_400_with_message(server):
    status, data = api(server, "/api/smiles", {"smiles": "C1CC"})
    assert status == 400 and not data["ok"] and "SMILES" in data["error"]
    status, data = api(server, "/api/optimize", {"molecule": {"atoms": "x"}})
    assert status == 400


def test_bad_json_and_unknown_route(server):
    port = server.server_address[1]
    conn = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
    conn.request("POST", "/api/info", body=b"{not json", headers={"Host": f"127.0.0.1:{port}", "X-Minimol": "1"})
    assert conn.getresponse().status == 400
    conn.close()
    status, _ = api(server, "/api/nope", {})
    assert status == 404
