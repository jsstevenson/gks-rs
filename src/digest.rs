use base64::{engine::general_purpose::URL_SAFE_NO_PAD, Engine as _};
use sha2::{Digest, Sha512};

pub fn sha512t24u(blob: &[u8]) -> String {
    let mut hasher = Sha512::new();
    hasher.update(blob);
    let result = hasher.finalize();
    let result_24 = &result[..24];
    URL_SAFE_NO_PAD.encode(result_24)
}
