/* np02_dxf.h  Newport 02 Library

Copyright (c) 2023 Mac Stevens <stevensm@earthlink.net> <www.macstevens.net>

Permission to use, copy, modify, and distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

Reference: https://opensource.org/licenses/ISC
*/
#ifndef NP02_DXF_H
#define NP02_DXF_H

#include <map>
#include <string>
#include <vector>
#include <stdint.h>

namespace np02 {



struct np02_dxf_file_init_params{
public:
    int32_t   init_val;         // TBD
};

enum np02_dxf_color{
    NP02_DXF_COLOR_000 = 0,
    NP02_DXF_COLOR_001_RED = 1,
    NP02_DXF_COLOR_002_YELLOW = 2,
    NP02_DXF_COLOR_003_GREEN = 3,
    NP02_DXF_COLOR_004_CYAN = 4,
    NP02_DXF_COLOR_005_BLUE = 5,
    NP02_DXF_COLOR_006_MAGENTA = 6,
    NP02_DXF_COLOR_007_BLACK = 7,
    NP02_DXF_COLOR_008 = 8,
    NP02_DXF_COLOR_009 = 9,
    /* ... */
    NP02_DXF_COLOR_255 = 255,
    NP02_DXF_COLOR_COUNT = 256
};


/* 
example: https://rlaanemets.com/hello-world-from-dxf/rectangle.dxf
*/
class np02_dxf_file{
private:
    typedef std::vector<std::string> str_vec;
    typedef str_vec::const_iterator str_vec_citr;
    typedef std::map<std::string, uint8_t> str_color_map;
    typedef str_color_map::iterator str_color_map_itr;
    typedef str_color_map::const_iterator str_color_map_citr;
private:
    str_vec m_dxf_shape_str_vec;
    str_color_map m_str_color_map;
public:
    np02_dxf_file();
    np02_dxf_file(const np02_dxf_file_init_params& init_params);
   ~np02_dxf_file();
    void init(const np02_dxf_file_init_params& init_params);
    void draw_line(const std::string& layer,const double& x0,const double& y0, 
        const double& x1, const double& y1, const uint8_t& color);
    void draw_box(const std::string& layer, const double& x0,const double& y0, 
        const double& x1, const double& y1, const uint8_t& color);
    void draw_circle(const std::string& layer, const double& x_ctr,
        const double& y_ctr, const double& radius, const uint8_t& color);
    void draw_arc(const std::string& layer, const double& x_ctr,
        const double& y_ctr, const double& radius,
        const double& start_angle_deg, const double& end_angle_deg,
        const uint8_t& color);
    void draw_text(const std::string& layer,const std::string& text,const double& x0,
        const double& y0, const double& height, const double& rot_deg,
        const uint8_t& color);
    int write_file(const char *file_name) const;
    int insert_in_file(const char *input_file_name,
        const char *output_file_name) const;
private:
    void add_layer( const std::string& layer, const uint8_t& color);
    void write_head( std::ostream& os )const;
    void write_shapes( std::ostream& os )const;
    void write_tail( std::ostream& os )const;
    static bool is_match_str(const std::string& x, const std::string& y);

};


class np02_dxf_test{
public:
    static int run_dxf_test();
};


}

#endif /* NP02_DXF_H */

