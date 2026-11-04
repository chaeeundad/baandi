fn main() {
    println!(
        "{}",
        serde_json::to_string_pretty(&baandi_mcp::schema::tools()).unwrap()
    );
}
