// image.hpp - minimal PNG writer and slice compositor.
//
// The Python renders with matplotlib. Nothing equivalent is linked here, and
// the ports must not depend on a Python toolchain, so figures are written
// directly as PNG. Pixel data is stored with deflate "stored" blocks, which
// needs no compression library: the files are a little larger than a normal
// PNG but are ordinary, readable PNGs.
#pragma once

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

namespace dhj {

inline std::uint32_t crc32_of(const std::uint8_t* data, std::size_t n, std::uint32_t crc = 0) {
    static std::uint32_t table[256];
    static bool ready = false;
    if (!ready) {
        for (std::uint32_t i = 0; i < 256; ++i) {
            std::uint32_t c = i;
            for (int k = 0; k < 8; ++k) c = (c & 1) ? 0xEDB88320u ^ (c >> 1) : c >> 1;
            table[i] = c;
        }
        ready = true;
    }
    crc = ~crc;
    for (std::size_t i = 0; i < n; ++i) crc = table[(crc ^ data[i]) & 0xFF] ^ (crc >> 8);
    return ~crc;
}

inline void put_u32_be(std::vector<std::uint8_t>& v, std::uint32_t x) {
    v.push_back(static_cast<std::uint8_t>(x >> 24));
    v.push_back(static_cast<std::uint8_t>(x >> 16));
    v.push_back(static_cast<std::uint8_t>(x >> 8));
    v.push_back(static_cast<std::uint8_t>(x));
}

inline void png_chunk(std::ofstream& out, const char tag[4], const std::vector<std::uint8_t>& data) {
    std::vector<std::uint8_t> hdr;
    put_u32_be(hdr, static_cast<std::uint32_t>(data.size()));
    out.write(reinterpret_cast<const char*>(hdr.data()), static_cast<std::streamsize>(hdr.size()));

    std::vector<std::uint8_t> body(tag, tag + 4);
    body.insert(body.end(), data.begin(), data.end());
    out.write(reinterpret_cast<const char*>(body.data()), static_cast<std::streamsize>(body.size()));

    std::vector<std::uint8_t> crc;
    put_u32_be(crc, crc32_of(body.data(), body.size()));
    out.write(reinterpret_cast<const char*>(crc.data()), static_cast<std::streamsize>(crc.size()));
}

// rgb must hold width * height * 3 bytes, row-major from the top.
inline bool write_png(const std::string& path, int width, int height,
                      const std::vector<std::uint8_t>& rgb) {
    std::ofstream out(path, std::ios::binary);
    if (!out) return false;

    static const std::uint8_t sig[8] = {0x89, 'P', 'N', 'G', 0x0D, 0x0A, 0x1A, 0x0A};
    out.write(reinterpret_cast<const char*>(sig), 8);

    std::vector<std::uint8_t> ihdr;
    put_u32_be(ihdr, static_cast<std::uint32_t>(width));
    put_u32_be(ihdr, static_cast<std::uint32_t>(height));
    ihdr.push_back(8);  // bit depth
    ihdr.push_back(2);  // colour type: truecolour RGB
    ihdr.push_back(0);  // deflate
    ihdr.push_back(0);  // adaptive filtering
    ihdr.push_back(0);  // no interlace
    png_chunk(out, "IHDR", ihdr);

    // Raw scanlines, each prefixed with filter type 0.
    std::vector<std::uint8_t> raw;
    raw.reserve(static_cast<std::size_t>(height) * (1 + static_cast<std::size_t>(width) * 3));
    for (int y = 0; y < height; ++y) {
        raw.push_back(0);
        const std::size_t off = static_cast<std::size_t>(y) * width * 3;
        raw.insert(raw.end(), rgb.begin() + static_cast<std::ptrdiff_t>(off),
                   rgb.begin() + static_cast<std::ptrdiff_t>(off + static_cast<std::size_t>(width) * 3));
    }

    // zlib stream: header, stored deflate blocks, adler32.
    std::vector<std::uint8_t> z;
    z.reserve(raw.size() + raw.size() / 65535 * 5 + 16);
    z.push_back(0x78);
    z.push_back(0x01);
    std::size_t pos = 0;
    while (pos < raw.size()) {
        const std::size_t len = std::min<std::size_t>(65535, raw.size() - pos);
        const bool final_block = (pos + len == raw.size());
        z.push_back(final_block ? 1 : 0);
        z.push_back(static_cast<std::uint8_t>(len & 0xFF));
        z.push_back(static_cast<std::uint8_t>(len >> 8));
        const std::uint16_t nlen = static_cast<std::uint16_t>(~len);
        z.push_back(static_cast<std::uint8_t>(nlen & 0xFF));
        z.push_back(static_cast<std::uint8_t>(nlen >> 8));
        z.insert(z.end(), raw.begin() + static_cast<std::ptrdiff_t>(pos),
                 raw.begin() + static_cast<std::ptrdiff_t>(pos + len));
        pos += len;
    }
    std::uint32_t a = 1, b = 0;
    for (std::uint8_t byte : raw) {
        a = (a + byte) % 65521;
        b = (b + a) % 65521;
    }
    put_u32_be(z, (b << 16) | a);
    png_chunk(out, "IDAT", z);

    png_chunk(out, "IEND", {});
    return static_cast<bool>(out);
}

// Simple RGB canvas with the few primitives the figures need.
struct Canvas {
    int w = 0, h = 0;
    std::vector<std::uint8_t> px;

    Canvas(int width, int height) : w(width), h(height), px(static_cast<std::size_t>(width) * height * 3, 255) {}

    void set(int x, int y, int r, int g, int b) {
        if (x < 0 || y < 0 || x >= w || y >= h) return;
        const std::size_t i = (static_cast<std::size_t>(y) * w + x) * 3;
        px[i] = static_cast<std::uint8_t>(r);
        px[i + 1] = static_cast<std::uint8_t>(g);
        px[i + 2] = static_cast<std::uint8_t>(b);
    }

    void rect(int x0, int y0, int x1, int y1, int r, int g, int b) {
        for (int y = y0; y < y1; ++y)
            for (int x = x0; x < x1; ++x) set(x, y, r, g, b);
    }

    // Circle outline, drawn thick enough to read at figure scale.
    void circle(double cx, double cy, double rx, double ry, int r, int g, int b, int thickness) {
        const int steps = 2048;
        for (int i = 0; i < steps; ++i) {
            const double t = 6.283185307179586 * i / steps;
            const double x = cx + rx * std::cos(t);
            const double y = cy + ry * std::sin(t);
            for (int dy = -thickness; dy <= thickness; ++dy)
                for (int dx = -thickness; dx <= thickness; ++dx)
                    set(static_cast<int>(x) + dx, static_cast<int>(y) + dy, r, g, b);
        }
    }

    // Dashed circle outline, for the target set.
    void dashed_circle(double cx, double cy, double rx, double ry, int r, int g, int b, int thickness) {
        const int steps = 2048;
        for (int i = 0; i < steps; ++i) {
            if ((i / 40) % 2) continue;
            const double t = 6.283185307179586 * i / steps;
            const double x = cx + rx * std::cos(t);
            const double y = cy + ry * std::sin(t);
            for (int dy = -thickness; dy <= thickness; ++dy)
                for (int dx = -thickness; dx <= thickness; ++dx)
                    set(static_cast<int>(x) + dx, static_cast<int>(y) + dy, r, g, b);
        }
    }
};

}  // namespace dhj
