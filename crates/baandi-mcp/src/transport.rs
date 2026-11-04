use crate::{Server, schema};
use serde_json::{Value, json};
use std::io::{self, BufRead, Read, Write};
const VERSION: &str = "2025-06-18";
fn rpc_error(id: Value, code: i32, message: &str) -> Value {
    json!({"jsonrpc":"2.0","id":id,"error":{"code":code,"message":message}})
}
impl Server {
    /// Process one UTF-8 JSON-RPC frame. Notifications never invoke tools.
    pub fn handle(&mut self, frame: &str) -> Option<Value> {
        if self.requests >= self.quotas.requests {
            return Some(rpc_error(Value::Null, -32000, "request quota exhausted"));
        }
        self.requests += 1;
        if frame.len() > self.quotas.frame_bytes {
            return Some(rpc_error(Value::Null, -32600, "frame quota exceeded"));
        }
        let v: Value = match serde_json::from_str(frame) {
            Ok(v) => v,
            Err(_) => return Some(rpc_error(Value::Null, -32700, "parse error")),
        };
        let id = v.get("id").cloned();
        let valid_id = id
            .as_ref()
            .is_none_or(|v| v.is_string() || v.as_i64().is_some() || v.as_u64().is_some());
        if !v.is_object()
            || v["jsonrpc"] != "2.0"
            || !v["method"].is_string()
            || !valid_id
            || v.get("params").is_some_and(|p| !p.is_object())
        {
            return Some(rpc_error(Value::Null, -32600, "invalid request"));
        }
        let method = v["method"].as_str().unwrap();
        if id.is_none() {
            if method == "notifications/initialized" && self.initialized {
                self.ready = true;
            }
            return None;
        }
        let id = id.unwrap();
        let result = match method {
            "initialize" => {
                if self.initialized {
                    return Some(rpc_error(id, -32600, "already initialized"));
                }
                if !v["params"]["protocolVersion"].is_string()
                    || !v["params"]["capabilities"].is_object()
                    || !v["params"]["clientInfo"]["name"].is_string()
                    || !v["params"]["clientInfo"]["version"].is_string()
                {
                    return Some(rpc_error(id, -32602, "invalid initialize parameters"));
                }
                self.initialized = true;
                json!({"protocolVersion":VERSION,"capabilities":{"tools":{"listChanged":false}},"serverInfo":{"name":"baandi-mcp","version":env!("CARGO_PKG_VERSION")},"instructions":"Observe or Supervised only. Host registers validated inputs and approves interventions outside MCP. No filesystem or shell tools."})
            }
            "ping" => json!({}),
            _ if !self.ready => {
                return Some(rpc_error(id, -32002, "initialize handshake required"));
            }
            "tools/list" => {
                if v["params"].get("cursor").is_some() {
                    return Some(rpc_error(id, -32602, "no pagination cursor supported"));
                }
                json!({"tools":schema::tools()})
            }
            "tools/call" => {
                let Some(name) = v["params"]["name"].as_str() else {
                    return Some(rpc_error(id, -32602, "tool name required"));
                };
                let specs = schema::tools();
                let Some(spec) = specs.as_array().unwrap().iter().find(|t| t["name"] == name)
                else {
                    return Some(rpc_error(id, -32602, "unknown tool"));
                };
                let args = v["params"]["arguments"].clone();
                if !schema::validate(&spec["inputSchema"], &args) {
                    return Some(rpc_error(id, -32602, "invalid tool arguments"));
                }
                match self.call(name, args) {
                    Ok(value) => {
                        json!({"content":[{"type":"text","text":value.to_string()}],"isError":false})
                    }
                    Err(message) => {
                        json!({"content":[{"type":"text","text":message}],"isError":true})
                    }
                }
            }
            _ => return Some(rpc_error(id, -32601, "method not found")),
        };
        Some(json!({"jsonrpc":"2.0","id":id,"result":result}))
    }
    /// Synchronous stdio framing with bounded input allocation. Oversized frames
    /// close the connection; never resume parsing the remainder as a new request.
    pub fn serve(&mut self, mut input: impl BufRead, mut output: impl Write) -> io::Result<()> {
        loop {
            if self.requests >= self.quotas.requests {
                return Ok(());
            }
            let mut frame = Vec::new();
            let n = (&mut input)
                .take(self.quotas.frame_bytes as u64 + 1)
                .read_until(b'\n', &mut frame)?;
            if n == 0 {
                return Ok(());
            }
            if n > self.quotas.frame_bytes {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "frame quota exceeded",
                ));
            }
            let reply = match std::str::from_utf8(&frame) {
                Ok(frame) => self.handle(frame),
                Err(_) => {
                    self.requests += 1;
                    Some(rpc_error(Value::Null, -32700, "invalid UTF-8"))
                }
            };
            if let Some(reply) = reply {
                serde_json::to_writer(&mut output, &reply)?;
                output.write_all(b"\n")?;
                output.flush()?;
            }
        }
    }
    pub fn serve_stdio(&mut self) -> io::Result<()> {
        self.serve(io::stdin().lock(), io::stdout().lock())
    }
}
