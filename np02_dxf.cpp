/* np02_dxf.cpp  Newport 02 Library

Copyright (c) 2023 Mac Stevens <stevensm@earthlink.net> <www.macstevens.net>

Permission to use, copy, modify, and distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

Reference: https://opensource.org/licenses/ISC
*/
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "np02_dxf.h"

namespace np02 {


np02_dxf_file::np02_dxf_file():m_dxf_shape_str_vec(){
}

np02_dxf_file::np02_dxf_file(const np02_dxf_file_init_params& init_params):
    m_dxf_shape_str_vec(){
init(init_params);
}

np02_dxf_file::~np02_dxf_file(){
}

void np02_dxf_file::init(const np02_dxf_file_init_params& init_params){
m_dxf_shape_str_vec.clear();
}

void np02_dxf_file::draw_line(const std::string& layer,const double& x0,
    const double& y0, const double& x1, const double& y1,const uint8_t& color){
std::ostringstream os;
os <<  "  0\n"
    << "LINE\n"
    << "  8\n"
    << layer << "\n"
    << " 10\n"
    << x0 << "\n"
    << " 20\n"
    << y0 << "\n"
    << " 11\n"
    << x1 << "\n"
    << " 21\n"
    << y1 << "\n"
    << " 62\n"
    << static_cast<int>(color) << "\n";
m_dxf_shape_str_vec.push_back(os.str());
add_layer(layer, color);
}

void np02_dxf_file::draw_box(const std::string& layer,const double& x0,
    const double& y0, const double& x1, const double& y1,const uint8_t& color){
draw_line(layer,x0,y0,x1,y0,color);
draw_line(layer,x1,y0,x1,y1,color);
draw_line(layer,x1,y1,x0,y1,color);
draw_line(layer,x0,y1,x0,y0,color);
}


void np02_dxf_file::draw_circle(const std::string& layer,
    const double& x_ctr,const double& y_ctr, 
    const double& radius, const uint8_t& color){
std::ostringstream os;
os <<  "  0\n"
    << "CIRCLE\n"
    << "  8\n"
    << layer << "\n"
    << " 10\n"
    << x_ctr << "\n"
    << " 20\n"
    << y_ctr << "\n"
    << " 40\n"
    << radius << "\n"
    << " 62\n"
    << static_cast<int>(color) << "\n";
m_dxf_shape_str_vec.push_back(os.str());
add_layer(layer, color);
}

void np02_dxf_file::draw_text(const std::string& layer,
    const std::string& text, const double& x0, const double& y0,
    const double& height, const double& rot_deg, const uint8_t& color){
std::ostringstream os;
os <<  "  0\n"
    << "TEXT\n"
    << "  1\n"
    << text << "\n"
    << "  8\n"
    << layer << "\n"
    << " 10\n"
    << x0 << "\n"
    << " 20\n"
    << y0 << "\n"
    << " 40\n"
    << height << "\n"
    << " 50\n"
    << rot_deg << "\n"
    << " 62\n"
    << static_cast<int>(color) << "\n";
m_dxf_shape_str_vec.push_back(os.str());
add_layer(layer, color);
}

int np02_dxf_file::write_file(const char *file_name) const{
int err_cnt = 0;
std::ofstream os(file_name, std::ios::out);

/* HEADER*/
os<<"  0\n"
    "SECTION\n"
    "  2\n"
    "HEADER\n"
    "  9\n"
    "$ACADVER\n"
    "  1\n"
    "AC1009\n"
    "  0\n"
    "ENDSEC\n";

#define NP02_DXF_PRINT_TABLES (0)
if(NP02_DXF_PRINT_TABLES){
    /* TABLES */
    os<<"  0\n"
        "SECTION\n"
        "  2\n"
        "TABLES\n"
        "  0\n"
        "TABLE\n"
        "  2\n"
        "LAYER\n";
    
    str_color_map_citr map_itr = m_str_color_map.begin();
    for(; map_itr != m_str_color_map.end(); ++map_itr){
        const std::string& lyr_name = map_itr->first;
        const uint8_t& color = map_itr->second;
        os<<"  0\n"
            "LAYER\n"
            "  2\n"
            << lyr_name << "\n"
            " 62\n"
            << (int)color << "\n"
            " 70\n"
            "     0\n";
        }
    
    os<<"  0\n" 
        "ENDTAB\n"
        "  0\n" 
        "ENDSEC\n";
    }

/* ENTITIES*/
os<<"  0\n"
    "SECTION\n"
    "  2\n"
    "ENTITIES\n";

str_vec_citr str_itr = m_dxf_shape_str_vec.begin();
for( ; str_itr != m_dxf_shape_str_vec.end(); ++str_itr){
    os << (*str_itr);
    }

os<<"  0\n" 
    "ENDSEC\n";

/* EOF */ 
os<<"  0\n" 
    "EOF\n";

if(!os.good()){ ++err_cnt; }

os.close();

return err_cnt;
}

void np02_dxf_file::add_layer( const std::string& layer,
    const uint8_t& color){
str_color_map_itr map_itr = m_str_color_map.lower_bound(layer);
if( ( map_itr == m_str_color_map.end() ) || ( map_itr->first != layer ) ){
    m_str_color_map.insert( map_itr, str_color_map::value_type(layer, color) );
    }
}


}
