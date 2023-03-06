/* np02_bmp.h  Newport 02 Library

Copyright (c) 2023 Mac Stevens <stevensm@earthlink.net> <www.macstevens.net>

Permission to use, copy, modify, and distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

Reference: https://opensource.org/licenses/ISC
*/
#ifndef NP02_BMP_H
#define NP02_BMP_H

#include <vector>
#include <stdint.h>

namespace np02 {

#if defined(__GNUC__) && (__GNUC__ < 5)
    #define PRAGMA_PACK_PUSH_POP_ENABLED (0)
#else
    #define PRAGMA_PACK_PUSH_POP_ENABLED (1)  
#endif

#if PRAGMA_PACK_PUSH_POP_ENABLED
    #pragma pack(push,2) 
#else
    #pragma pack(2)
#endif
struct np02_bmp_header{         // Total: 54 bytes
public:
    uint16_t  type;             // Magic identifier: 0x4d42
    uint32_t  size;             // File size in bytes
    uint16_t  reserved1;        // Not used
    uint16_t  reserved2;        // Not used
    uint32_t  offset;           // Offset to image data in bytes from beginning of file (54 bytes)
    uint32_t  dib_header_size;  // DIB Header size in bytes (40 bytes)
    int32_t   width_px;         // Width of the image
    int32_t   height_px;        // Height of image
    uint16_t  num_planes;       // Number of color planes
    uint16_t  bits_per_pixel;   // Bits per pixel
    uint32_t  compression;      // Compression type
    uint32_t  image_size_bytes; // Image size in bytes
    int32_t   x_resolution_ppm; // Pixels per meter
    int32_t   y_resolution_ppm; // Pixels per meter
    uint32_t  num_colors;       // Number of colors  
    uint32_t  important_colors; // Important colors 
};
#if PRAGMA_PACK_PUSH_POP_ENABLED
    #pragma pack(pop) 
#else
    #pragma pack(4)
#endif

struct np02_bmp_file_init_params{
public:
    int32_t   width_px;         // Width of the image
    int32_t   height_px;        // Height of image
};

struct np02_bmp_color{
public:
    typedef np02_bmp_color color_t;
public: /* https://en.wikipedia.org/wiki/Web_colors#Basic_colors */
    static np02_bmp_color white  (){return color_t(0xFF,0xFF,0xFF);}
    static np02_bmp_color silver (){return color_t(0xC0,0xC0,0xC0);}
    static np02_bmp_color gray   (){return color_t(0x80,0x80,0x80);}
    static np02_bmp_color black  (){return color_t(0x00,0x00,0x00);}
    static np02_bmp_color red    (){return color_t(0xFF,0x00,0x00);}
    static np02_bmp_color maroon (){return color_t(0x80,0x00,0x00);}
    static np02_bmp_color yellow (){return color_t(0xFF,0xFF,0x00);}
    static np02_bmp_color olive  (){return color_t(0x80,0x80,0x00);}
    static np02_bmp_color lime   (){return color_t(0x00,0xFF,0x00);}
    static np02_bmp_color green  (){return color_t(0x00,0x80,0x00);}
    static np02_bmp_color aqua   (){return color_t(0x00,0xFF,0xFF);}
    static np02_bmp_color teal   (){return color_t(0x00,0x80,0x80);}
    static np02_bmp_color blue   (){return color_t(0x00,0x00,0xFF);}
    static np02_bmp_color navy   (){return color_t(0x00,0x00,0x80);}
    static np02_bmp_color fuchsia(){return color_t(0xFF,0x00,0xFF);}
    static np02_bmp_color purple (){return color_t(0x80,0x00,0x80);}

    static color_t bl_white  (){return color_t(0xFF,0xFF,0xFF,true);}
    static color_t bl_silver (){return color_t(0xC0,0xC0,0xC0,true);}
    static color_t bl_gray   (){return color_t(0x80,0x80,0x80,true);}
    static color_t bl_black  (){return color_t(0x00,0x00,0x00,true);}
    static color_t bl_red    (){return color_t(0xFF,0x00,0x00,true);}
    static color_t bl_maroon (){return color_t(0x80,0x00,0x00,true);}
    static color_t bl_yellow (){return color_t(0xFF,0xFF,0x00,true);}
    static color_t bl_olive  (){return color_t(0x80,0x80,0x00,true);}
    static color_t bl_lime   (){return color_t(0x00,0xFF,0x00,true);}
    static color_t bl_green  (){return color_t(0x00,0x80,0x00,true);}
    static color_t bl_aqua   (){return color_t(0x00,0xFF,0xFF,true);}
    static color_t bl_teal   (){return color_t(0x00,0x80,0x80,true);}
    static color_t bl_blue   (){return color_t(0x00,0x00,0xFF,true);}
    static color_t bl_navy   (){return color_t(0x00,0x00,0x80,true);}
    static color_t bl_fuchsia(){return color_t(0xFF,0x00,0xFF,true);}
    static color_t bl_purple (){return color_t(0x80,0x00,0x80,true);}
public:
    uint8_t m_red;
    uint8_t m_green;
    uint8_t m_blue;
    bool    m_blend;
public:
    np02_bmp_color(uint8_t r, uint8_t g, uint8_t b, bool blend=false):
        m_red(r),m_green(g),m_blue(b),m_blend(blend){}
};


class np02_bmp_file{
private:
    typedef std::vector<uint8_t> uint8_vec;
    typedef uint8_vec::iterator uint8_vec_itr;
    typedef uint8_vec::const_iterator uint8_vec_citr;
    typedef std::vector<uint8_vec> uint8_vec_vec;
    typedef uint8_vec_vec::iterator uint8_vec_vec_itr;
    typedef uint8_vec_vec::const_iterator uint8_vec_vec_citr;
private:
    static const int m_hershey_simplex_font_coord[95][112];
    np02_bmp_header m_header;
    uint8_vec_vec m_rows;
public:
    np02_bmp_file();
    np02_bmp_file(const np02_bmp_file_init_params& init_params);
   ~np02_bmp_file();
    void init(const np02_bmp_file_init_params& init_params);
    void draw_pixel(const int32_t& i,const int32_t& j, const np02_bmp_color& color);
    void draw_line(const int32_t& i0,const int32_t& j0, 
        const int32_t& i1,const int32_t& j1, const np02_bmp_color& color);
    void draw_box(const int32_t& i0,const int32_t& j0, 
        const int32_t& i1,const int32_t& j1, const np02_bmp_color& color);
    void draw_diamond(const int32_t& i_ctr,const int32_t& j_ctr, 
        const int32_t& w, const np02_bmp_color& color);
    void draw_wide_line(const double& ii0, const double& jj0, 
        const double& ii1,const double& jj1, const double& ww,
        const np02_bmp_color& color);
    void draw_circle(const double& ii_ctr, const double& jj_ctr, 
        const double& rr, const np02_bmp_color& color);
    void draw_text(const char *text, const double& ii0, const double& jj0, 
        const double& hh, const double& rot_deg, const np02_bmp_color& color);
    void write_file(const char *file_name) const;
};


}

#endif /* NP02_BMP_H */

