import { readFileSync } from "fs";
import { dirname, join } from "path";
import { fileURLToPath } from "url";
import init from "../pkg/moyo_wasm.js";

const __dirname = dirname(fileURLToPath(import.meta.url));
const wasm_binary_path = join(__dirname, "..", "pkg", "moyo_wasm_bg.wasm");

// Initialize WASM module for all tests
const wasm_bytes = readFileSync(wasm_binary_path);
await init(wasm_bytes);
