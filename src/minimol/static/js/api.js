// Calls to the local minimol server. Every call is a POST with a JSON body and
// the X-Minimol header the server requires.

export async function api(path, body = {}) {
  let res;
  try {
    res = await fetch(path, {
      method: "POST",
      headers: { "Content-Type": "application/json", "X-Minimol": "1" },
      body: JSON.stringify(body),
    });
  } catch {
    throw new Error("Cannot reach the minimol server. Is it still running?");
  }
  let data;
  try {
    data = await res.json();
  } catch {
    throw new Error(`Server error (HTTP ${res.status})`);
  }
  if (!data.ok) throw new Error(data.error || `Server error (HTTP ${res.status})`);
  return data;
}
