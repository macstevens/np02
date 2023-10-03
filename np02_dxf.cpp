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
    << "100\n"
    << "AcDbEntity\n"
    << "  8\n"
    << layer << "\n"
    << " 62\n"
    << static_cast<int>(color) << "\n"
    << "100\n"
    << "AcDbLine\n"
    << " 10\n"
    << x0 << "\n"
    << " 20\n"
    << y0 << "\n"
    << " 11\n"
    << x1 << "\n"
    << " 21\n"
    << y1 << "\n";
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
    << "100\n"
    << "AcDbEntity\n"
    << "  8\n"
    << layer << "\n"
    << " 62\n"
    << static_cast<int>(color) << "\n"
    << "100\n"
    << "AcDbCircle\n"
    << " 10\n"
    << x_ctr << "\n"
    << " 20\n"
    << y_ctr << "\n"
    << " 40\n"
    << radius << "\n";
m_dxf_shape_str_vec.push_back(os.str());
add_layer(layer, color);
}

void np02_dxf_file::draw_arc(const std::string& layer, const double& x_ctr,
    const double& y_ctr, const double& radius,
    const double& start_angle_deg, const double& end_angle_deg,
    const uint8_t& color){
std::ostringstream os;
double start_angle = fmod( start_angle_deg, 360.0 );
if( start_angle < 0 ){ start_angle += 360.0; }
double end_angle = fmod( end_angle_deg, 360.0 );
if( end_angle < 0 ){ end_angle += 360.0; }
os <<  "  0\n"
    << "ARC\n"
    << "100\n"
    << "AcDbEntity\n"
    << "  8\n"
    << layer << "\n"
    << " 62\n"
    << static_cast<int>(color) << "\n"
    << "100\n"
    << "AcDbCircle\n"
    << " 10\n"
    << x_ctr << "\n"
    << " 20\n"
    << y_ctr << "\n"
    << " 40\n"
    << radius << "\n"
    << "100\n"
    << "AcDbArc\n"
    << " 50\n"
    << start_angle << "\n"
    << " 51\n"
    << end_angle << "\n";
m_dxf_shape_str_vec.push_back(os.str());
add_layer(layer, color);
}

void np02_dxf_file::draw_text(const std::string& layer,
    const std::string& text, const double& x0, const double& y0,
    const double& height, const double& rot_deg, const uint8_t& color){
std::ostringstream os;
os <<  "  0\n"
    << "TEXT\n"
    << "100\n"
    << "AcDbEntity\n"
    << "  8\n"
    << layer << "\n"
    << " 62\n"
    << static_cast<int>(color) << "\n"
    << "100\n"
    << "AcDbText\n"
    << "  1\n"
    << text << "\n"
    << " 10\n"
    << x0 << "\n"
    << " 20\n"
    << y0 << "\n"
    << " 40\n"
    << height << "\n"
    << " 50\n"
    << rot_deg << "\n";
m_dxf_shape_str_vec.push_back(os.str());
add_layer(layer, color);
}

int np02_dxf_file::write_file(const char *file_name) const{
int err_cnt = 0;
std::ofstream os(file_name, std::ios::out);
write_head(os);
write_shapes(os);
write_tail(os);
if(!os.good()){ ++err_cnt; }
os.close();
return err_cnt;
}

int np02_dxf_file::insert_in_file(const char *input_file_name,
    const char *output_file_name) const{
int err_cnt = 0;

/* TODO: if input file name == output file name, rename input file
   before accessing and remove the renamed file afterwards. */

/* open input & output files */
std::ifstream ifs(input_file_name);
if(!ifs.good()){ ++err_cnt; }
std::ofstream ofs(output_file_name, std::ios::out);
if(!ofs.good()){ ++err_cnt; }

if( 0 == err_cnt ){
    std::string line;
    const size_t rcnt_hist_sz = 4;
    std::vector<std::string> recent_lines(rcnt_hist_sz);
    size_t rcnt_hist_idx = 0;
    const std::string *line1 = NULL;
    const std::string *line2 = NULL;
    const std::string *line3 = NULL;
    const std::string *line4 = NULL;

    /* copy lines up to and including beginning of ENTITIES section */
    bool entities_header_found = false;
    while( ifs.good() && ofs.good() && !entities_header_found ){
        /* read line */
        getline(ifs,line);

        /* copy to output file*/
        ofs << line << "\n";

        /* update recent history */
        recent_lines.at(rcnt_hist_idx) = line;
        line1 = line2;
        line2 = line3;
        line3 = line4;
        line4 = &(recent_lines.at(rcnt_hist_idx));
        rcnt_hist_idx = (rcnt_hist_idx+1) % rcnt_hist_sz;

        /* check for ENTITIES section beginning */
        if(NULL != line1){
            assert(!entities_header_found);
            assert(NULL != line2);
            assert(NULL != line3);
            assert(NULL != line4);
            static const std::string hdr_line1 = "  0";
            static const std::string hdr_line2 = "SECTION";
            static const std::string hdr_line3 = "  2";
            static const std::string hdr_line4 = "ENTITIES";
            if( is_match_str(*line1, hdr_line1) &&
                is_match_str(*line2, hdr_line2) &&
                is_match_str(*line3, hdr_line3) &&
                is_match_str(*line4, hdr_line4) ){
                entities_header_found = true;
                }
            }
        }

    /* insert shapes before pre-existing shapes */
    write_shapes(ofs);

    /* copy remainder of file */
    while( ifs.good() && ofs.good() ){
        /* read line */
        getline(ifs,line);

        /* copy to output file*/
        ofs << line << "\n";
        }

    if( !ofs.good() ){
        ++err_cnt;
        }
    }

return err_cnt;
}


void np02_dxf_file::add_layer( const std::string& layer,
    const uint8_t& color){
str_color_map_itr map_itr = m_str_color_map.lower_bound(layer);
if( ( map_itr == m_str_color_map.end() ) || ( map_itr->first != layer ) ){
    m_str_color_map.insert( map_itr, str_color_map::value_type(layer, color) );
    }
}

void np02_dxf_file::write_head( std::ostream& os ) const{
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
}

void np02_dxf_file::write_shapes( std::ostream& os ) const{
str_vec_citr str_itr = m_dxf_shape_str_vec.begin();
for( ; str_itr != m_dxf_shape_str_vec.end(); ++str_itr){
    os << (*str_itr);
    }
}

void np02_dxf_file::write_tail( std::ostream& os ) const{
os<<"  0\n" 
    "ENDSEC\n";

/* EOF */ 
os<<"  0\n" 
    "EOF\n";
}

/* compare case-insensitive until final non-whitespace */
bool np02_dxf_file::is_match_str(const std::string& x, const std::string& y){
static const std::string whitespace = " \t\r\n\v\f";
bool match = true;
const std::string::size_type last_pos_x = x.find_last_not_of(whitespace);
const std::string::size_type last_pos_y = y.find_last_not_of(whitespace);
if( (last_pos_x == std::string::npos ) ||
    (last_pos_y == std::string::npos ) ||
    (last_pos_x != last_pos_y ) ){
    match = false;
    }
else{
    std::string::const_iterator x_itr = x.begin(); 
    const std::string::const_iterator x_end_itr = x_itr + (last_pos_x + 1);
    std::string::const_iterator y_itr = y.begin(); 
    for(; match && (x_itr != x_end_itr); ++x_itr, ++y_itr ){
        assert( x_itr != x.end() );
        assert( y_itr != y.end() );
        const char x_ch = toupper(*x_itr);
        const char y_ch = toupper(*y_itr);
        if( x_ch != y_ch ){
            match = false;
            }
        }
    }
return match;
}


int np02_dxf_test::run_dxf_test(){
int err_cnt = 0;

/* first shape set */
np02_dxf_file_init_params p0;
p0.init_val = 0;
np02_dxf_file dxf_file_0(p0);
dxf_file_0.draw_line("A", 0.0, 0.0, 3.0, 4.0, NP02_DXF_COLOR_001_RED );
dxf_file_0.draw_box("A", 0.0, 0.0, 3.0, 4.0, NP02_DXF_COLOR_002_YELLOW );
dxf_file_0.draw_circle("B", 0.0, 0.0, 5.0, NP02_DXF_COLOR_003_GREEN );
dxf_file_0.draw_arc("B", 6.0, 8.0, 5.0, 135.0, 270.0, NP02_DXF_COLOR_004_CYAN);
dxf_file_0.draw_text("B", "TEXT ABC", -3.0, 4.0, 0.5, 0.0,
    NP02_DXF_COLOR_005_BLUE);

/* write to file */
static const std::string file_name_0 = "test_file_0.dxf";
err_cnt += dxf_file_0.write_file(file_name_0.c_str());

/* first shape set */
np02_dxf_file_init_params p1;
p1.init_val = 0;
np02_dxf_file dxf_file_1(p1);
dxf_file_1.draw_line("C", 5.0, 0.0, 8.0, 4.0, NP02_DXF_COLOR_005_BLUE );
dxf_file_1.draw_box("C", 5.0, 0.0, 8.0, 4.0, NP02_DXF_COLOR_005_BLUE );
dxf_file_1.draw_circle("D", 5.0, 0.0, 5.0, NP02_DXF_COLOR_006_MAGENTA );
dxf_file_1.draw_arc("D", 11.0, 8.0, 5.0, 135.0, 270.0,
    NP02_DXF_COLOR_006_MAGENTA);
dxf_file_1.draw_text("D", "TEXT XYZ", 2.0, 4.0, 0.5, 0.0,
    NP02_DXF_COLOR_006_MAGENTA);

/* write to file */
static const std::string file_name_1 = "test_file_1.dxf";
err_cnt += dxf_file_1.write_file(file_name_1.c_str());

/* write a third file with both contents */
static const std::string file_name_12 = "test_file_12.dxf";
err_cnt += dxf_file_1.insert_in_file( file_name_0.c_str(),
    file_name_12.c_str() );

/* clean up */
err_cnt += remove(file_name_0.c_str());
err_cnt += remove(file_name_1.c_str());
err_cnt += remove(file_name_12.c_str());

return err_cnt;
}


} /* namespace np02 */
