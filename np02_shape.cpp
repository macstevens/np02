/* np02_shape.cpp  Newport 02 Library

Copyright (c) 2023 Mac Stevens <stevensm@earthlink.net> <www.macstevens.net>

Permission to use, copy, modify, and distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

Reference: https://opensource.org/licenses/ISC
*/

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <limits>
#include <set>
#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "np02.h"
#include "np02_bmp.h"


#define CF01_SUPPORT ( 1 )
#if defined( CF01_SUPPORT )
#include "cf01.h"
    #define AA_INCR_CALL_DEPTH()        CF01_AA_INCR_CALL_DEPTH()
    #define AA_DECR_CALL_DEPTH()        CF01_AA_DECR_CALL_DEPTH()
    #define AUTO_ASSERT( _condition )   CF01_AUTO_ASSERT(_condition)
    #define AA_ALWAYS_ASSERT( _condition ) \
                  CF01_AA_XDBG_ASSERT((_condition), CF01_AA_DEBUG_LEVEL_0)
    #define AA_SHOULD_RUN_XDBG(_dbg_lvl) CF01_AA_SHOULD_RUN_XDBG(_dbg_lvl)
    #define AA_XDBG_ASSERT( _condition, _dbg_lvl ) \
                            CF01_AA_XDBG_ASSERT( (_condition), _dbg_lvl )
    #define AA_ERR_BUF()                CF01_AA_ERR_BUF()
    #define AA_ERR_BUF_CAPACITY()       CF01_AA_ERR_BUF_CAPACITY()
    #define AA_ERR_BUF_POS_PTR()        CF01_AA_ERR_BUF_POS_PTR()
    #define AA_DEBUG_LEVEL_0            CF01_AA_DEBUG_LEVEL_0
    #define AA_DEBUG_LEVEL_1            CF01_AA_DEBUG_LEVEL_1
    #define AA_DEBUG_LEVEL_2            CF01_AA_DEBUG_LEVEL_2
    #define AA_DEBUG_LEVEL_3            CF01_AA_DEBUG_LEVEL_3
#else
    #define AA_INCR_CALL_DEPTH()        
    #define AA_DECR_CALL_DEPTH()        
    #define AUTO_ASSERT( _condition )   assert(_condition)
    #define AA_ALWAYS_ASSERT( _condition )   assert(_condition)
    #define AA_SHOULD_RUN_XDBG(_dbg_lvl)  (0)
    #define AA_XDBG_ASSERT( _condition, _dbg_lvl ) 
    #define AA_ERR_BUF()                (NULL)
    #define AA_ERR_BUF_CAPACITY()       (0)
    #define AA_ERR_BUF_POS_PTR()        (NULL)
    #define AA_DEBUG_LEVEL_0            (0)
    #define AA_DEBUG_LEVEL_1            (1)
    #define AA_DEBUG_LEVEL_2            (2)
    #define AA_DEBUG_LEVEL_3            (3)
#endif

namespace np02 {


np02_shape::np02_shape():m_shape_type(NP02_SHAPE_TYPE_SHAPE),
    m_shape_idx(0), m_shp_owner_idx(NP02_SHP_OWNER_INVALID_IDX),
    m_shp_alloc(NULL),  m_head_loc_grid_node(NULL){}

/* free resources */
void np02_shape::destruct(){
if(NULL != m_head_loc_grid_node ){
    np02_loc_grid *loc_grid = m_head_loc_grid_node->get_loc_grid();
    /* m_head_loc_grid_node is used as free chain pointer, so
    m_head_loc_grid_node->get_loc_grid() might be NULL if 
    m_shp_owner_idx is invalid */
    AA_ALWAYS_ASSERT( (NULL != loc_grid) ||
       (NP02_SHP_OWNER_INVALID_IDX == m_shp_owner_idx) )
    if(NULL != loc_grid){
        loc_grid->remove_shape_from_loc_grid(this);
        }
    }
}

np02_shape::~np02_shape(){
destruct();
}

lyr_idx_type np02_shape::get_lyr_idx() const{
return (NULL == m_head_loc_grid_node) ? NP02_LYR_INVALID_IDX : 
    m_head_loc_grid_node->get_lyr_idx();
}

np02_loc_grid *np02_shape::get_loc_grid() const{
np02_loc_grid *g = (NULL == m_head_loc_grid_node) ? NULL :
    m_head_loc_grid_node->get_loc_grid();
return g;
}

void np02_shape::get_local_shapes(type_idx_shp_vec *local_shapes) const{
AA_ALWAYS_ASSERT(local_shapes != NULL);
const np02_loc_grid_node *n = get_head_loc_grid_node();
for( ; NULL != n; n = n->get_s_next() ){
    const np02_loc_grid_node *m = n->get_prev();
    for( ; NULL != m; m = m->get_prev() ){
        np02_shape *s = m->get_owner();
        AA_ALWAYS_ASSERT(s != NULL);
        AA_ALWAYS_ASSERT(s != this);
        const np02_shape::type_idx_pair sort_key(
            s->get_shape_type(), s->get_shape_idx());
        const np02_shape::type_idx_shp sort_entry(sort_key, s);
        local_shapes->push_back(sort_entry);   
        }
    m = n->get_next();
    for( ; NULL != m; m = m->get_next() ){
        np02_shape *s = m->get_owner();
        AA_ALWAYS_ASSERT(s != NULL);
        AA_ALWAYS_ASSERT(s != this);
        const np02_shape::type_idx_pair sort_key(
            s->get_shape_type(), s->get_shape_idx());
        const np02_shape::type_idx_shp sort_entry(sort_key, s);
        local_shapes->push_back(sort_entry);   
        }
    }
std::sort(local_shapes->begin(), local_shapes->end());
type_idx_shp_vec_itr end_pos = std::unique(local_shapes->begin(),
    local_shapes->end());
local_shapes->erase( end_pos, local_shapes->end() );
}

double np02_shape::get_distance_from_shape_double_check(
    const np02_shape *s, double *err_estimate, np02_xy *near_xy,
    np02_xy *other_near_xy) const{
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != s);

/* get overall bounding box */
np02_xy bb_xy_min;
np02_xy bb_xy_max;
get_bb(&bb_xy_min, &bb_xy_max);
np02_xy s_bb_xy_min;
np02_xy s_bb_xy_max;
s->get_bb(&s_bb_xy_min, &s_bb_xy_max);
np02_xy xy_min = bb_xy_min;
np02_xy xy_max = bb_xy_max;
if(s_bb_xy_min.get_x() < xy_min.get_x()){
    xy_min.set_x(s_bb_xy_min.get_x()); }
if(s_bb_xy_min.get_y() < xy_min.get_y()){
    xy_min.set_y(s_bb_xy_min.get_y()); }
if(s_bb_xy_max.get_x() > xy_max.get_x()){
    xy_max.set_x(s_bb_xy_max.get_x()); }
if(s_bb_xy_max.get_y() > xy_max.get_y()){
    xy_max.set_y(s_bb_xy_max.get_y()); }

np02_xy search_ctr( (xy_min.get_x() + xy_max.get_x())/2.0,
                        (xy_min.get_y() + xy_max.get_y())/2.0 );
double search_half_range_x = xy_max.get_x() - search_ctr.get_x();
double search_half_range_y = xy_max.get_y() - search_ctr.get_y();
double d_near = std::numeric_limits<double>::max();
np02_xy xy_near = search_ctr;
np02_xy this_nr_xy, other_nr_xy;
np02_xy xy(0.0, 0.0);
for(size_t h = 0; h < 5; ++h) {
    const double x_step = search_half_range_x / 2.0;
    const double y_step = search_half_range_y / 2.0;
    for(double i = -2.0; i <= 2.0; i += 1.0){
        xy.set_x(search_ctr.get_x() + (i * x_step));
        for(double j = -2.0; j <= 2.0; j += 1.0){
            xy.set_y(search_ctr.get_y() + (i * y_step));
            np02_xy xy_near_a, xy_near_b;
            const double d_a = get_distance_from_xy(xy, &xy_near_a);
            const double d_b = s->get_distance_from_xy(xy, &xy_near_b);
            const double d = d_a + d_b;
            if(d < d_near){
                d_near = d;
                xy_near = xy;
                this_nr_xy = xy_near_a;
                other_nr_xy = xy_near_b;
                }
            }
        }
    search_half_range_x /= 4.0;
    search_half_range_y /= 4.0;
    }
d_near = this_nr_xy.get_distance_to(other_nr_xy);
if(NULL != err_estimate){
    *err_estimate = 2.0 * (search_half_range_x + search_half_range_y);
    }
if(NULL != near_xy){ *near_xy = this_nr_xy; }
if(NULL != other_near_xy){ *other_near_xy = other_nr_xy; }

AA_DECR_CALL_DEPTH();
return d_near;
}

void np02_shape::translate(const np02_xy& dxy){
AA_INCR_CALL_DEPTH();
np02_loc_grid *loc_grid = get_loc_grid();
if(NULL != loc_grid){ loc_grid->remove_shape_from_loc_grid(this); }
AUTO_ASSERT(NULL == get_loc_grid());
translate_no_loc_grid(dxy);
if(NULL != loc_grid){ loc_grid->insert_shape_in_loc_grid(this); }
AUTO_ASSERT(get_loc_grid() == loc_grid);
AA_DECR_CALL_DEPTH();
}

void np02_shape::rotate(const np02_xy& rot_ctr, const double& rot_deg){
AA_INCR_CALL_DEPTH();
np02_loc_grid *loc_grid = get_loc_grid();
if(NULL != loc_grid){ loc_grid->remove_shape_from_loc_grid(this); }
AUTO_ASSERT(NULL == get_loc_grid());
rotate_no_loc_grid(rot_ctr, rot_deg);
if(NULL != loc_grid){ loc_grid->insert_shape_in_loc_grid(this); }
AUTO_ASSERT(get_loc_grid() == loc_grid);
AA_DECR_CALL_DEPTH();
}


int np02_shape::verify_distance_from_shape_result(
    const np02_shape *shape_a,
    const np02_shape *shape_b, const np02_xy& xy_near_a,
    const np02_xy& xy_near_b, const double& distance_from,
    char *err_msg, const size_t err_msg_capacity, size_t *err_msg_pos ){
int err_cnt = 0;

if(NULL == shape_a ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }
if(NULL == shape_b ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }

if((NULL != shape_a) && (NULL != shape_b)){
    /* bounding box */
    const double xy_near_a_b_d = xy_near_a.get_distance_to(xy_near_b);
    const double xy_near_a_b_d_err=fabs(xy_near_a_b_d-fabs(distance_from));

    /* bounding box */
    np02_xy bb_a_xy_min, bb_a_xy_max;
    np02_xy bb_b_xy_min, bb_b_xy_max;
    np02_xy bb_ab_xy_min, bb_ab_xy_max;
    shape_a->get_bb(&bb_a_xy_min, &bb_a_xy_max);
    shape_b->get_bb(&bb_b_xy_min, &bb_b_xy_max);
    bb_ab_xy_min.set_x( ( bb_a_xy_min.get_x() < bb_b_xy_min.get_x() ) ?
                          bb_a_xy_min.get_x() : bb_b_xy_min.get_x() );
    bb_ab_xy_min.set_y( ( bb_a_xy_min.get_y() < bb_b_xy_min.get_y() ) ?
                          bb_a_xy_min.get_y() : bb_b_xy_min.get_y() );
    bb_ab_xy_max.set_x( ( bb_a_xy_max.get_x() > bb_b_xy_max.get_x() ) ?
                          bb_a_xy_max.get_x() : bb_b_xy_max.get_x() );
    bb_ab_xy_max.set_y( ( bb_a_xy_max.get_y() > bb_b_xy_max.get_y() ) ?
                          bb_a_xy_max.get_y() : bb_b_xy_max.get_y() );
    
    /* bounding box overlap 
           AAAAAAAAAA                             AAAAAAAAAA        
           min    max                             min    max        
                                                     
                   BBBBBBBBBB              BBBBBBBBBB
                   min    max              min    max
    */
    const bool bb_overlap = ( ( bb_a_xy_min.get_x() <= bb_b_xy_max.get_x() ) &&
                              ( bb_b_xy_min.get_x() <= bb_a_xy_max.get_x() ) &&
                              ( bb_a_xy_min.get_y() <= bb_b_xy_max.get_y() ) &&
                              ( bb_b_xy_min.get_y() <= bb_a_xy_max.get_y() ) );
    
    /* max err distance */
    const double bb_ab_semi_perimeter = 
        ( bb_ab_xy_max.get_x() - bb_ab_xy_min.get_x() ) +
        ( bb_ab_xy_max.get_y() - bb_ab_xy_min.get_y() );
    const double max_err_d = bb_ab_semi_perimeter * 1e-6;

    /* distances to xy_near_a, xy_near_b 

      +------------+                                   +------------+
      |            |                                   |            |
      |    A       |                                   |            |
      |            |   d>0    +--------------+         |     A      |
      |   xy_near_a*----------*xy_near_b     |         |        +--------+
      |            |          |              |         |        |d<0|    |
      +------------+          |      B       |         +--------|---+    |
                              |              |                  |    B   |
                              |              |                  +--------+
                              |              |
                              +--------------+
    */
    np02_xy xy_near_a_a(0.0,0.0);
    const double d_a_a = shape_a->get_distance_from_xy(xy_near_a,&xy_near_a_a);
    np02_xy xy_near_a_b(0.0,0.0);
    const double d_a_b = shape_a->get_distance_from_xy(xy_near_b,&xy_near_a_b);
    np02_xy xy_near_b_a(0.0,0.0);
    const double d_b_a = shape_b->get_distance_from_xy(xy_near_a,&xy_near_b_a);
    np02_xy xy_near_b_b(0.0,0.0);
    const double d_b_b = shape_b->get_distance_from_xy(xy_near_b,&xy_near_b_b);

    /* checks */
    if((distance_from > max_err_d) && ( xy_near_a_b_d_err > max_err_d) ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s  a_type=%i  b_type=%i"
            "  distance_from=%f > %f  &&  xy_near_a_b_d_err=%f > %f", 
            __FILE__, __LINE__, __FUNCTION__,
            shape_a->get_shape_type(), shape_b->get_shape_type(),
            distance_from, max_err_d, xy_near_a_b_d_err, max_err_d );
        }

    if( (!bb_overlap) && (distance_from < -max_err_d) ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s  a_type=%i  b_type=%i"
            "  (!bb_overlap) && distance_from=%f < %f", 
            __FILE__, __LINE__, __FUNCTION__,
            shape_a->get_shape_type(), shape_b->get_shape_type(),
            distance_from, -max_err_d );
        }

    /*   */
    if( d_a_a > max_err_d ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s  a_type=%i  b_type=%i"
            "  d_a_a=%f > %f", 
            __FILE__, __LINE__, __FUNCTION__,
            shape_a->get_shape_type(), shape_b->get_shape_type(),
            d_a_a, max_err_d );
        }
    /* if( (distance_from < -max_err_d) && (d_a_b > max_err_d) ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s  a_type=%i  b_type=%i"
            "  distance_from=%f < %f  &&  d_a_b=%f > %f", 
            __FILE__, __LINE__, __FUNCTION__,
            shape_a->get_shape_type(), shape_b->get_shape_type(),
            distance_from, -max_err_d, d_a_b, max_err_d );
        } */
    if( (distance_from > max_err_d) && (d_a_b < -max_err_d) ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s  a_type=%i  b_type=%i"
            "  distance_from=%f > %f  &&  d_a_b=%f < %f", 
            __FILE__, __LINE__, __FUNCTION__,
            shape_a->get_shape_type(), shape_b->get_shape_type(),
            distance_from, max_err_d, d_a_b, -max_err_d );
        }

    if( d_b_b > max_err_d ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s  a_type=%i  b_type=%i"
            "  d_b_b=%f > %f", 
            __FILE__, __LINE__, __FUNCTION__,
            shape_a->get_shape_type(), shape_b->get_shape_type(),
            d_b_b, max_err_d );
        }
    /* if( (distance_from < -max_err_d) && (d_b_a > max_err_d) ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s  a_type=%i  b_type=%i"
            "  distance_from=%f < %f  &&  d_b_a=%f > %f", 
            __FILE__, __LINE__, __FUNCTION__,
            shape_a->get_shape_type(), shape_b->get_shape_type(),
            distance_from, -max_err_d, d_b_a, max_err_d );
        } */
    if( (distance_from > max_err_d) && (d_b_a < -max_err_d) ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s  a_type=%i  b_type=%i"
            "  distance_from=%f > %f  &&  d_b_a=%f < %f", 
            __FILE__, __LINE__, __FUNCTION__,
            shape_a->get_shape_type(), shape_b->get_shape_type(),
            distance_from, max_err_d, d_b_a, -max_err_d );
        }
    }

return err_cnt;
}


int np02_shape::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const{
int err_cnt = 0;
if((m_shape_type < NP02_SHAPE_TYPE_SHAPE ) ||
    ( m_shape_type > NP02_SHAPE_TYPE_COUNT)){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shape: this=%x  m_shape_type=%i\n", this, m_shape_type); 
    }
if(m_shape_idx > NP02_SHAPE_MAX_IDX){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shape: this=%x  m_shape_idx=%i\n", this, m_shape_idx); 
    }

/* m_head_loc_grid_node is used as free chain pointer, so
m_head_loc_grid_node->get_loc_grid() might be NULL if 
m_shp_owner_idx is invalid */
if(( NULL != m_head_loc_grid_node) && 
    (NP02_SHP_OWNER_INVALID_IDX == m_shp_owner_idx) ){
    const np02_loc_grid *loc_grid = m_head_loc_grid_node->get_loc_grid();
    np02_loc_grid_dim loc_grid_dim;
    if(NULL == loc_grid){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shape: this=%x head_loc_grid_node(=%x)->loc_grid == NULL\n",
            this, m_head_loc_grid_node);
        loc_grid_dim.reset();
        }
    else{
        loc_grid_dim = loc_grid->get_loc_grid_dim();
        }
    size_t max_loc_grid_node_count = 
        static_cast<size_t>(loc_grid_dim.get_w()) *
        static_cast<size_t>(loc_grid_dim.get_h());
    if(0 == max_loc_grid_node_count){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shape: this=%x head_loc_grid_node(=%x)->loc_grid(=%x)"
            "->loc_grid_dim  w=%i  h=%i\n",
            this, m_head_loc_grid_node, loc_grid,loc_grid_dim.get_w(),
            loc_grid_dim.get_h() );
        }
    size_t loc_grid_node_count = 1;
    const np02_loc_grid_node *lg_node = m_head_loc_grid_node;
    while((NULL!=lg_node) && (loc_grid_node_count<=max_loc_grid_node_count)){
        err_cnt += lg_node->verify_data(err_msg,err_msg_capacity,err_msg_pos);
        if(this != lg_node->get_owner()){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shape: this=%x  !=  (lg_node=%x)->owner=%x\n",
                this, lg_node, lg_node->get_owner() );
            }
        if(loc_grid != lg_node->get_loc_grid()){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shape: this=%x  loc_grid=%x  !=  (lg_node=%x)->loc_grid=%x\n",
                this, loc_grid, lg_node, lg_node->get_loc_grid() );
            }
        const np02_loc_grid_node *s_next_lg_node = lg_node->get_s_next();
        if(NULL != s_next_lg_node){
            ++loc_grid_node_count;
            if( lg_node != s_next_lg_node->get_s_prev()){
                ++err_cnt;
                np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                    "shape: this=%x lg_node=%x != (s_next_lg_node=%x)"
                    "->s_prev=%x\n", this, lg_node, s_next_lg_node,
                    s_next_lg_node->get_s_prev() );
                }
            }
        lg_node = s_next_lg_node;
        }
    if(loc_grid_node_count > max_loc_grid_node_count){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shape: this=%x  loc_grid_node_count=%i > "
            "max_loc_grid_node_count=%i\n",
            this, loc_grid_node_count, max_loc_grid_node_count );
        }
    }
return err_cnt; 
}

std::ostream& np02_shape::ostream_output(std::ostream& os) const{
os << "<shape>\n";
os << std::hex;
os << "<this>" << this << "</this>\n";
os << std::dec;
os << "<shape_type>" << static_cast<int>(m_shape_type) << "</shape_type>\n";
os << "<shape_idx>" << m_shape_idx << "</shape_idx>\n";
os << "<shp_owner_idx>" << m_shp_owner_idx << "</shp_owner_idx>\n";
os << std::hex;
os << "<shp_alloc>" << m_shp_alloc << "</shp_alloc>\n";
os << "<head_loc_grid_node>" << m_head_loc_grid_node << "</head_loc_grid_node>\n";
os << std::dec;
os << "</shape>\n";
return os;
}

void np02_shape::write_bmp_file(const np02_xy& xy_min,
    const double& pixel_num, const np02_bmp_color& color,
    np02_bmp_file *bmp_file) const{
np02_xy bb_xy_min, bb_xy_max;
get_bb(&bb_xy_min, &bb_xy_max);
const int32_t i_min = static_cast<int32_t>(
    (bb_xy_min.get_x() - xy_min.get_x()) * pixel_num);
const int32_t i_max = static_cast<int32_t>(
    (bb_xy_max.get_x() - xy_min.get_x()) * pixel_num);
const int32_t j_min = static_cast<int32_t>(
    (bb_xy_min.get_y() - xy_min.get_y()) * pixel_num);
const int32_t j_max = static_cast<int32_t>(
    (bb_xy_max.get_y() - xy_min.get_y()) * pixel_num);

np02_xy xy;
for(int32_t i = i_min; i <= i_max; ++i){
    xy.set_x(xy_min.get_x() + (static_cast<double>(i) / pixel_num));
    for(int32_t j = j_min; j <= j_max; ++j){
        xy.set_y(xy_min.get_y() + (static_cast<double>(j) / pixel_num));
        const double d = get_distance_from_xy(xy);
        if(d <= 0.0 ){
            bmp_file->draw_pixel(i,j,color);
            }
        }
    }
}

void np02_shape::write_dxf_file(const std::string& layer,
    const uint8_t& color, np02_dxf_file *dxf_file) const{

}


np02_circle::np02_circle():m_ctr(0.0, 0.0), m_radius(0.0){
set_shape_type(NP02_SHAPE_TYPE_CIRCLE);
}

np02_circle::~np02_circle(){}

void np02_circle::get_bb(np02_xy *xy_min, np02_xy *xy_max) const{
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != xy_min);
AA_ALWAYS_ASSERT(NULL != xy_max);
xy_min->set_x(m_ctr.get_x()-m_radius);
xy_min->set_y(m_ctr.get_y()-m_radius);
xy_max->set_x(m_ctr.get_x()+m_radius);
xy_max->set_y(m_ctr.get_y()+m_radius);
AA_DECR_CALL_DEPTH();
}

void np02_circle::get_loc_grid_indices_for_init(
    const np02_loc_grid_dim& loc_grid_dim, const double& extra_search_d,
    np02_uint16_pair_vec *index_vec) const{
AA_INCR_CALL_DEPTH();
AUTO_ASSERT(m_radius >= 0.0);
AUTO_ASSERT(loc_grid_dim.get_sq_size() > 0.0);
AUTO_ASSERT(extra_search_d >= 0.0);
AA_ALWAYS_ASSERT(NULL != index_vec);
uint16_t i,j;
np02_xy xy_min, xy_max;
get_bb(&xy_min, &xy_max);
xy_min.set_x(xy_min.get_x()-extra_search_d);
xy_min.set_y(xy_min.get_y()-extra_search_d);
xy_max.set_x(xy_max.get_x()+extra_search_d);
xy_max.set_y(xy_max.get_y()+extra_search_d);
np02_uint16_pair ij_min, ij_max;
loc_grid_dim.get_bb_indices(xy_min, xy_max, &ij_min, &ij_max);
const double d_threshold = m_radius + extra_search_d + 
    (0.75 * loc_grid_dim.get_sq_size());
AUTO_ASSERT(d_threshold > 0.0);
const double dsq_threshold = d_threshold * d_threshold;
bool i_on_loc_grid_edge, j_on_loc_grid_edge, should_add;
for( i = ij_min.first; i <= ij_max.first; ++i ){
    const double x = loc_grid_dim.get_sq_ctr_x(i);
    i_on_loc_grid_edge = ((0 == i) || (loc_grid_dim.get_w() == (i+1))) ? 
        true : false;
    for( j = ij_min.second; j <= ij_max.second; ++j ){
        j_on_loc_grid_edge = ((0 == j) || (loc_grid_dim.get_h() == (j+1))) ?
            true : false;
        if(i_on_loc_grid_edge || j_on_loc_grid_edge){
            should_add = true;
            }
        else{
            const double y = loc_grid_dim.get_sq_ctr_y(j);
            const double dx = x - m_ctr.get_x();
            const double dy = y - m_ctr.get_y();
            const double dsq = (dx*dx) + (dy*dy);
            should_add = (dsq <= dsq_threshold) ? true : false;
            }
        if(should_add){
            index_vec->push_back(np02_uint16_pair(i,j));
            }
        }
    }
AA_DECR_CALL_DEPTH();
}

double np02_circle::get_distance_from_xy(const np02_xy& xy,
    np02_xy *near_xy) const{
const double dx = xy.get_x() - m_ctr.get_x();
const double dy = xy.get_y() - m_ctr.get_y();
const double dsq = (dx*dx) + (dy*dy);
const double d_to_ctr = sqrt(dsq);
const double d = d_to_ctr - m_radius;
if(NULL != near_xy){
    if(m_radius < 1e-40) {
        *near_xy = m_ctr; }
    else if(d_to_ctr < 1e-40) {
        if(fabs(dx) > fabs(dy)){
            near_xy->set_x(m_ctr.get_x() + ((dx > 0) ? m_radius : -m_radius));
            near_xy->set_y(0.0);
            }
        else{
            near_xy->set_x(0.0);
            near_xy->set_y(m_ctr.get_y() + ((dy > 0) ? m_radius : -m_radius));
            }
        }
    else{ 
        const double f = m_radius / d_to_ctr;
        near_xy->set_x(m_ctr.get_x() + (f * dx));
        near_xy->set_y(m_ctr.get_y() + (f * dy));
        }
    }
return d;
}

double np02_circle::get_distance_from_line_seg_ab(const np02_xy& xy_a,
    const np02_xy& xy_b, np02_xy *near_xy,
    np02_xy *other_near_xy) const{
AA_INCR_CALL_DEPTH();
AUTO_ASSERT(m_radius >= 0.0);
double dx, dy, dsq, d;
double d_ctr = 0.0;

/* forward unit vector A->B */
const double dx_ab = xy_b.get_x() - xy_a.get_x();
const double dy_ab = xy_b.get_y() - xy_a.get_y();
const double dsq_ab = (dx_ab * dx_ab) + (dy_ab * dy_ab);
np02_xy fwd_ab(1.0,0.0);
if(dsq_ab > 1e-50){
    const double d_ab = sqrt(dsq_ab);
    fwd_ab.set_x(dx_ab/d_ab);
    fwd_ab.set_y(dy_ab/d_ab);
    }
else if(fabs(dx_ab) > fabs(dy_ab)){
    fwd_ab.set_x(dx_ab > 0.0 ? 1.0 : -1.0);
    fwd_ab.set_y(0.0);
    }
else{
    fwd_ab.set_x(0.0);
    fwd_ab.set_y(dy_ab > 0.0 ? 1.0 : -1.0);
    }

const double fwd_dot_ctr = (fwd_ab.get_x() * m_ctr.get_x()) + 
                           (fwd_ab.get_y() * m_ctr.get_y());
const double fwd_dot_a = (fwd_ab.get_x() * xy_a.get_x()) + 
                         (fwd_ab.get_y() * xy_a.get_y());
const double fwd_dot_b = (fwd_ab.get_x() * xy_b.get_x()) + 
                         (fwd_ab.get_y() * xy_b.get_y());

if(fwd_dot_ctr<fwd_dot_a){
    /* Point A is the closest point to circle center */
    dx = xy_a.get_x() - m_ctr.get_x();
    dy = xy_a.get_y() - m_ctr.get_y();
    dsq = (dx*dx) + (dy*dy);
    d_ctr = sqrt(dsq);
    if(NULL != near_xy){
        if(m_radius < 1e-40) {
            *near_xy = m_ctr; }
        else if(d_ctr < 1e-40) {
            if(fabs(dx) > fabs(dy)){
                near_xy->set_x(m_ctr.get_x() + ((dx > 0) ? m_radius : -m_radius));
                near_xy->set_y(0.0);
                }
            else{
                near_xy->set_x(0.0);
                near_xy->set_y(m_ctr.get_y() + ((dy > 0) ? m_radius : -m_radius));
                }
            }
        else{ 
            const double f = m_radius / d_ctr;
            near_xy->set_x(m_ctr.get_x() + (f * dx));
            near_xy->set_y(m_ctr.get_y() + (f * dy));
            }
        }
    if(NULL != other_near_xy){*other_near_xy = xy_a; }
    }
else if(fwd_dot_ctr>fwd_dot_b){
    /* Point B is the closest point to circle center */
    dx = xy_b.get_x() - m_ctr.get_x();
    dy = xy_b.get_y() - m_ctr.get_y();
    dsq = (dx*dx) + (dy*dy);
    d_ctr = sqrt(dsq);
    if(NULL != near_xy){
        if(m_radius < 1e-40) {
            *near_xy = m_ctr; }
        else if(d_ctr < 1e-40) {
            if(fabs(dx) > fabs(dy)){
                near_xy->set_x(m_ctr.get_x() + ((dx > 0) ? m_radius : -m_radius));
                near_xy->set_y(0.0);
                }
            else{
                near_xy->set_x(0.0);
                near_xy->set_y(m_ctr.get_y() + ((dy > 0) ? m_radius : -m_radius));
                }
            }
        else{ 
            const double f = m_radius / d_ctr;
            near_xy->set_x(m_ctr.get_x() + (f * dx));
            near_xy->set_y(m_ctr.get_y() + (f * dy));
            }
        }
    if(NULL != other_near_xy){*other_near_xy = xy_b; }
    }
else{
    /* circle center is closest to a point on line segment between A and B. */
    const np02_xy xy_ab((xy_a.get_x() + xy_b.get_x())/2.0,
                            (xy_a.get_y() + xy_b.get_y())/2.0);
    const double fwd_cross_ab = (fwd_ab.get_x() * xy_ab.get_y()) - 
                                (fwd_ab.get_y() * xy_ab.get_x());
    const double fwd_cross_ctr = (fwd_ab.get_x() * m_ctr.get_y()) - 
                                 (fwd_ab.get_y() * m_ctr.get_x());
    d_ctr = fabs(fwd_cross_ab - fwd_cross_ctr);

    if(NULL != near_xy){
        const double fwd_cross_circle_near_xy = fwd_cross_ctr - m_radius;
        near_xy->set_x( ( fwd_dot_ctr * fwd_ab.get_x() ) - 
                        ( fwd_cross_circle_near_xy * fwd_ab.get_y() ) );
        near_xy->set_y( ( fwd_dot_ctr * fwd_ab.get_y() ) + 
                        ( fwd_cross_circle_near_xy * fwd_ab.get_x() ) );
        }
    if(NULL != other_near_xy){
        other_near_xy->set_x( ( fwd_dot_ctr * fwd_ab.get_x() ) - 
                              ( fwd_cross_ab * fwd_ab.get_y() ) );
        other_near_xy->set_y( ( fwd_dot_ctr * fwd_ab.get_y() ) + 
                              ( fwd_cross_ab * fwd_ab.get_x() ) );
        }
    }
d = d_ctr - m_radius;
AA_DECR_CALL_DEPTH();
return d;
}

double np02_circle::get_distance_from_circle(const np02_circle *c,
    np02_xy *near_xy, np02_xy *circle_near_xy) const{
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != c);
AUTO_ASSERT(m_radius >= 0.0);
AUTO_ASSERT(c->get_radius() >= 0.0);
const np02_xy& c_ctr = c->get_ctr();
const double dx = c_ctr.get_x() - m_ctr.get_x();
const double dy = c_ctr.get_y() - m_ctr.get_y();
const double dsq = (dx*dx) + (dy*dy);
const double d_ctr = sqrt(dsq);
const double d = d_ctr - (m_radius + c->get_radius());
if((NULL != near_xy) || (NULL != circle_near_xy)){
    np02_xy unit_v(1.0,0.0);
    if( d_ctr < 1e-40){
        if(fabs(dx) > fabs(dy)){
            unit_v.set_x((dx>0.0) ? 1.0:-1.0);
            unit_v.set_y(0.0);
            }
        else{
            unit_v.set_x(0.0);
            unit_v.set_y((dy>0.0) ? 1.0:-1.0);
            }
        }
    else{
        unit_v.set_x(dx/d_ctr);
        unit_v.set_y(dy/d_ctr);
        }
    if(NULL != near_xy){
        near_xy->set_x(m_ctr.get_x() + (m_radius * unit_v.get_x()));
        near_xy->set_y(m_ctr.get_y() + (m_radius * unit_v.get_y()));
        }
    if(NULL != circle_near_xy){
        circle_near_xy->set_x(c_ctr.get_x() -
            ((c->get_radius()) * unit_v.get_x()));
        circle_near_xy->set_y(c_ctr.get_y() -
            ((c->get_radius()) * unit_v.get_y()));
        }
    }
AA_DECR_CALL_DEPTH();
return d;
}

double np02_circle::get_distance_from_line_seg(
    const np02_line_seg *n, np02_xy *near_xy,
    np02_xy *line_seg_near_xy) const{
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != n);
AUTO_ASSERT(m_radius >= 0.0);
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_1);
np02_xy seg_near_xy;
const double d_ctr = n->get_distance_from_xy(m_ctr, &seg_near_xy);
const double d = d_ctr - m_radius;
if(NULL != near_xy){
    if( d_ctr < 0.0 ){
        *near_xy = m_ctr;
        }
    else if (d < 0.0 ){
        *near_xy = seg_near_xy;
        }
    else {
        const double dx = seg_near_xy.get_x() - m_ctr.get_x();
        const double dy = seg_near_xy.get_y() - m_ctr.get_y();
        const double ddsq = (dx*dx) + (dy*dy);
        const double dd = sqrt(ddsq);
        AUTO_ASSERT(fabs(dd - fabs(d_ctr)) <= fabs((dd + fabs(d_ctr))*1e-6) );
        np02_xy unit_v(1.0,0.0);
        if( dd < 1e-40){
            if(fabs(dx) > fabs(dy)){
                unit_v.set_x((dx>0.0) ? 1.0:-1.0);
                unit_v.set_y(0.0);
                }
            else{
                unit_v.set_x(0.0);
                unit_v.set_y((dy>0.0) ? 1.0:-1.0);
                }
            }
        else{
            unit_v.set_x(dx/dd);
            unit_v.set_y(dy/dd);
            }
        near_xy->set_x(m_ctr.get_x() + (m_radius * unit_v.get_x()));
        near_xy->set_y(m_ctr.get_y() + (m_radius * unit_v.get_y()));
        }
    }
if(NULL != line_seg_near_xy){*line_seg_near_xy = seg_near_xy; }
AA_DECR_CALL_DEPTH();
return d;
}

double np02_circle::get_distance_from_rect(const np02_rect *r,
    np02_xy *near_xy, np02_xy *rect_near_xy) const{
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != r);
AUTO_ASSERT(m_radius >= 0.0);
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_1);
np02_xy rectangle_near_xy;
const double d_ctr = r->get_distance_from_xy(m_ctr, &rectangle_near_xy);
const double d = d_ctr - m_radius;
if(NULL != near_xy){
    if( d > 0.0 ){
        const double dx = rectangle_near_xy.get_x() - m_ctr.get_x();
        const double dy = rectangle_near_xy.get_y() - m_ctr.get_y();
        const double ddsq = (dx*dx) + (dy*dy);
        const double dd = sqrt(ddsq);
        AUTO_ASSERT(fabs(dd - fabs(d_ctr)) <= fabs((dd + fabs(d_ctr))*1e-6) );
        np02_xy unit_v(1.0,0.0);
        if( dd < 1e-40){
            if(fabs(dx) > fabs(dy)){
                unit_v.set_x((dx>0.0) ? 1.0:-1.0);
                unit_v.set_y(0.0);
                }
            else{
                unit_v.set_x(0.0);
                unit_v.set_y((dy>0.0) ? 1.0:-1.0);
                }
            }
        else{
            unit_v.set_x(dx/dd);
            unit_v.set_y(dy/dd);
            }
        near_xy->set_x(m_ctr.get_x() + (m_radius * unit_v.get_x()));
        near_xy->set_y(m_ctr.get_y() + (m_radius * unit_v.get_y()));
        }
    else if( d_ctr > 0.0 ){
        *near_xy = rectangle_near_xy;
        }
    else{
        *near_xy = m_ctr;
        }
    }
if(NULL != rect_near_xy){*rect_near_xy = rectangle_near_xy; }
AA_DECR_CALL_DEPTH();
return d;
}

double np02_circle::get_distance_from_polygon(const np02_polygon *p,
    np02_xy *near_xy, np02_xy *polygon_near_xy) const{
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != p);
AUTO_ASSERT(m_radius >= 0.0);
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_1);
np02_xy poly_near_xy;
const double d_ctr = p->get_distance_from_xy(m_ctr, &poly_near_xy);
const double d = d_ctr - m_radius;
if(NULL != near_xy){
    if(d_ctr > 0.0){
        const double dx = poly_near_xy.get_x() - m_ctr.get_x();
        const double dy = poly_near_xy.get_y() - m_ctr.get_y();
        const double ddsq = (dx*dx) + (dy*dy);
        const double dd = sqrt(ddsq);
        AUTO_ASSERT(fabs(dd - d_ctr) <= fabs((dd + d_ctr)*1e-6) );
        np02_xy unit_v(1.0,0.0);
        if( dd < 1e-40){
            if(fabs(dx) > fabs(dy)){
                unit_v.set_x((dx>0.0) ? 1.0:-1.0);
                unit_v.set_y(0.0);
                }
            else{
                unit_v.set_x(0.0);
                unit_v.set_y((dy>0.0) ? 1.0:-1.0);
                }
            }
        else{
            unit_v.set_x(dx/dd);
            unit_v.set_y(dy/dd);
            }
        near_xy->set_x(m_ctr.get_x() + (m_radius * unit_v.get_x()));
        near_xy->set_y(m_ctr.get_y() + (m_radius * unit_v.get_y()));
        }
    else{
        *near_xy = m_ctr;
        }
    }
if(NULL != polygon_near_xy){*polygon_near_xy = poly_near_xy; }
AA_DECR_CALL_DEPTH();
return d;
}

double np02_circle::get_distance_from_shape(const np02_shape *s,
    np02_xy *near_xy, np02_xy *other_near_xy) const{
double d = 0.0;
AA_INCR_CALL_DEPTH();
if(NULL != s){
    const np02_circle *c = dynamic_cast<const np02_circle *>(s);
    const np02_line_seg *n = dynamic_cast<const np02_line_seg *>(s);
    const np02_rect *r = dynamic_cast<const np02_rect *>(s);
    const np02_polygon *p = dynamic_cast<const np02_polygon *>(s);
    if(NULL != c){ d = get_distance_from_circle(c, near_xy, other_near_xy); }
    else if(NULL!=n){ d=get_distance_from_line_seg(n, near_xy, other_near_xy);}
    else if(NULL!=r){ d=get_distance_from_rect(r, near_xy, other_near_xy); }
    else if(NULL!=p){ d=get_distance_from_polygon(p, near_xy, other_near_xy); }
    else{ AA_ALWAYS_ASSERT(false); }
    }
AA_DECR_CALL_DEPTH();
return d;
}

void np02_circle::translate_no_loc_grid(const np02_xy& dxy){
m_ctr.set_x( m_ctr.get_x() + dxy.get_x() );
m_ctr.set_y( m_ctr.get_y() + dxy.get_y() );
}

void np02_circle::rotate_no_loc_grid(const np02_xy& rot_ctr,
    const double& rot_deg){
const np02_xy rot_arm_initial(m_ctr.get_x() - rot_ctr.get_x(),
    m_ctr.get_y() - rot_ctr.get_y());
if( (rot_arm_initial.get_x() != 0.0) || (rot_arm_initial.get_y() != 0.0) ){
    double cos_rot = 1.0;
    double sin_rot = 0.0;
    if(0.0 == rot_deg){ cos_rot = 1.0; sin_rot = 0.0; }
    else if(90.0 == rot_deg){ cos_rot = 0.0; sin_rot = 1.0; }
    else if(-90.0 == rot_deg){ cos_rot = 0.0; sin_rot = -1.0; }
    else if(180.0 == rot_deg){ cos_rot = -1.0; sin_rot = 0.0; }
    else if(-180.0 == rot_deg){ cos_rot = -1.0; sin_rot = 0.0; }
    else{
        const double rot_rad =
            rot_deg * (3.1415926535897932384626433832795 / 180.0);
        cos_rot = cos(rot_rad);
        sin_rot = sin(rot_rad);
        }
    const np02_xy rot_arm_final(
        (rot_arm_initial.get_x() * cos_rot) -
        (rot_arm_initial.get_y() * sin_rot),
        (rot_arm_initial.get_x() * sin_rot) +
        (rot_arm_initial.get_y() * cos_rot));
    m_ctr.set_x(rot_ctr.get_x() + rot_arm_final.get_x());
    m_ctr.set_y(rot_ctr.get_y() + rot_arm_final.get_y());
    }
}

int np02_circle::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const{
int err_cnt = 0;
err_cnt += np02_shape::verify_data(err_msg, err_msg_capacity, err_msg_pos );

if(NP02_SHAPE_TYPE_CIRCLE != get_shape_type() ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "circle: this=%x  m_shape_type=%i\n", this, get_shape_type()); 
    }

const np02_shp_alloc *shp_alloc = get_shp_alloc();
if((NULL != shp_alloc) && 
    (this != shp_alloc->alloc_get_circle_by_idx(get_shape_idx()))){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "circle: this=%x  != (shp_alloc=%x) circle_by_idx(idx=%i)=%x\n",
        this, shp_alloc, get_shape_idx(),
        shp_alloc->alloc_get_circle_by_idx(get_shape_idx())); 
    }

if(m_radius < 0.0){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "circle: this=%x m_radius=%f < 0.0\n", this, m_radius);
    }
return err_cnt;
}

std::ostream& np02_circle::ostream_output(std::ostream& os) const{
os << "<circle>\n";
os << std::hex;
os << "<this>" << this << "</this>\n";
os << std::dec;
np02_shape::ostream_output(os);
os << "<ctr><x>" << m_ctr.get_x() << "</x>";
os << "<y>" << m_ctr.get_y() << "</y></ctr>\n";
os << "<radius>" << m_radius << "</radius>\n";
os << "</circle>\n";
return os;
}

void np02_circle::write_bmp_file(const np02_xy& xy_min,
    const double& pixel_num, const np02_bmp_color& color,
    np02_bmp_file *bmp_file) const{
np02_shape::write_bmp_file(xy_min, pixel_num, color, bmp_file);
}

void np02_circle::write_dxf_file(const std::string& layer,
    const uint8_t& color, np02_dxf_file *dxf_file) const{
AA_INCR_CALL_DEPTH();
AUTO_ASSERT(NULL != dxf_file);
dxf_file->draw_circle( layer, m_ctr.get_x(), m_ctr.get_y(), m_radius, color );
AA_DECR_CALL_DEPTH();
}


np02_rect::np02_rect():m_ctr(0.0,0.0),m_w(0.0),m_h(0.0),m_rot_deg(0.0),
    m_fwd(1.0,0.0), m_fwd_dot_ctr(0.0), m_fwd_cross_ctr(0.0),
    m_p00(0.0,0.0),m_p01(0.0,0.0),m_p10(0.0,0.0),m_p11(0.0,0.0){
set_shape_type(NP02_SHAPE_TYPE_RECT);
}

np02_rect::~np02_rect(){}

/* does not update loctor grid*/
void np02_rect::init(const np02_xy& ctr, const double& w, const double& h,
    const double& rot_deg){
AA_INCR_CALL_DEPTH();
AUTO_ASSERT(w >= 0.0);
AUTO_ASSERT(h >= 0.0);
m_ctr = ctr;
m_w = w;
m_h = h;
m_rot_deg = fmod(rot_deg, 360.0);

if(m_rot_deg == 0.0){ m_fwd.set_x(1.0); m_fwd.set_y(0.0); }
else if(m_rot_deg == 90.0){ m_fwd.set_x(0.0); m_fwd.set_y(1.0); }
else if(m_rot_deg == 180.0){ m_fwd.set_x(-1.0); m_fwd.set_y(0.0); }
else if(m_rot_deg == 270.0){ m_fwd.set_x(0.0); m_fwd.set_y(-1.0); }
else if(m_rot_deg == -90.0){ m_fwd.set_x(0.0); m_fwd.set_y(-1.0); }
else if(m_rot_deg == -180.0){ m_fwd.set_x(-1.0); m_fwd.set_y(0.0); }
else if(m_rot_deg == -270.0){ m_fwd.set_x(0.0); m_fwd.set_y(1.0); }
else{
    const double rot_rad = m_rot_deg*(3.1415926535897932384626433832795/180.0);
    m_fwd.set_x(cos(rot_rad));
    m_fwd.set_y(sin(rot_rad));
    }

/* p01                                     p11
    *---- ---- ---- ----^---- ---- ---- ----*
    |                   |                   |
    |                   |offset_cross       |
    |                   |                   |
    |                   |     offset_fwd    |
                        +------------------->                                   
    |                  ctr                  |       cross
    |                                       |       ^
    |                                       |       |
    |                                       |       |
    *---- ---- ---- ---- ---- ---- ---- ----*       +------>fwd
   p00                                     p10
*/

m_fwd_dot_ctr = (m_fwd.get_x()*m_ctr.get_x()) + (m_fwd.get_y()*m_ctr.get_y());
m_fwd_cross_ctr=(m_fwd.get_x()*m_ctr.get_y()) - (m_fwd.get_y()*m_ctr.get_x());
const np02_xy cross(-m_fwd.get_y(), m_fwd.get_x());
const np02_xy offset_fwd(m_w * m_fwd.get_x()/2.0, m_w * m_fwd.get_y()/2.0);
const np02_xy offset_cross(m_h*cross.get_x()/2.0, m_h * cross.get_y()/2.0);

const np02_xy offset_11( offset_fwd.get_x() + offset_cross.get_x(),
                             offset_fwd.get_y() + offset_cross.get_y() );
const np02_xy offset_10( offset_fwd.get_x() - offset_cross.get_x(),
                             offset_fwd.get_y() - offset_cross.get_y() );

m_p00.set_x( m_ctr.get_x() - offset_11.get_x() );
m_p00.set_y( m_ctr.get_y() - offset_11.get_y() );
m_p01.set_x( m_ctr.get_x() - offset_10.get_x() );
m_p01.set_y( m_ctr.get_y() - offset_10.get_y() );
m_p10.set_x( m_ctr.get_x() + offset_10.get_x() );
m_p10.set_y( m_ctr.get_y() + offset_10.get_y() );
m_p11.set_x( m_ctr.get_x() + offset_11.get_x() );
m_p11.set_y( m_ctr.get_y() + offset_11.get_y() );

AA_XDBG_ASSERT(0 == verify_data_num(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_1);
AA_DECR_CALL_DEPTH();
}

void np02_rect::get_bb(np02_xy *xy_min, np02_xy *xy_max) const{
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != xy_min);
AA_ALWAYS_ASSERT(NULL != xy_max);
AA_XDBG_ASSERT(0 == verify_data_num(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_2);
*xy_min = m_ctr;
*xy_max = m_ctr;
if( m_p00.get_x() < xy_min->get_x() ){ xy_min->set_x(m_p00.get_x()); }
if( m_p01.get_x() < xy_min->get_x() ){ xy_min->set_x(m_p01.get_x()); }
if( m_p10.get_x() < xy_min->get_x() ){ xy_min->set_x(m_p10.get_x()); }
if( m_p11.get_x() < xy_min->get_x() ){ xy_min->set_x(m_p11.get_x()); }
if( m_p00.get_y() < xy_min->get_y() ){ xy_min->set_y(m_p00.get_y()); }
if( m_p01.get_y() < xy_min->get_y() ){ xy_min->set_y(m_p01.get_y()); }
if( m_p10.get_y() < xy_min->get_y() ){ xy_min->set_y(m_p10.get_y()); }
if( m_p11.get_y() < xy_min->get_y() ){ xy_min->set_y(m_p11.get_y()); }
if( m_p00.get_x() > xy_max->get_x() ){ xy_max->set_x(m_p00.get_x()); }
if( m_p01.get_x() > xy_max->get_x() ){ xy_max->set_x(m_p01.get_x()); }
if( m_p10.get_x() > xy_max->get_x() ){ xy_max->set_x(m_p10.get_x()); }
if( m_p11.get_x() > xy_max->get_x() ){ xy_max->set_x(m_p11.get_x()); }
if( m_p00.get_y() > xy_max->get_y() ){ xy_max->set_y(m_p00.get_y()); }
if( m_p01.get_y() > xy_max->get_y() ){ xy_max->set_y(m_p01.get_y()); }
if( m_p10.get_y() > xy_max->get_y() ){ xy_max->set_y(m_p10.get_y()); }
if( m_p11.get_y() > xy_max->get_y() ){ xy_max->set_y(m_p11.get_y()); }
AA_DECR_CALL_DEPTH();
}

void np02_rect::get_loc_grid_indices_for_init(
    const np02_loc_grid_dim& loc_grid_dim,
    const double& extra_search_d, np02_uint16_pair_vec *index_vec)
    const{
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != index_vec);
AUTO_ASSERT(extra_search_d >= 0.0);
AA_XDBG_ASSERT(0 == loc_grid_dim.verify_data(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_2);
AA_XDBG_ASSERT(0 == verify_data_num(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_2);

uint16_t i,j;
np02_xy xy_min, xy_max;
get_bb(&xy_min, &xy_max);
xy_min.set_x(xy_min.get_x() - extra_search_d);
xy_min.set_y(xy_min.get_y() - extra_search_d);
xy_max.set_x(xy_max.get_x() + extra_search_d);
xy_max.set_y(xy_max.get_y() + extra_search_d);

np02_uint16_pair ij_min, ij_max;
loc_grid_dim.get_bb_indices(xy_min, xy_max, &ij_min, &ij_max);
const double d_threshold = extra_search_d + (0.75 * loc_grid_dim.get_sq_size());
AUTO_ASSERT(d_threshold > 0.0);

const double min_fwd_dot = get_fwd_dot_0() - d_threshold;
const double max_fwd_dot = get_fwd_dot_1() + d_threshold;
const double min_fwd_cross = get_fwd_cross_0() - d_threshold;
const double max_fwd_cross = get_fwd_cross_1() + d_threshold;

bool i_on_loc_grid_edge, j_on_loc_grid_edge, should_add;
for( i = ij_min.first; i <= ij_max.first; ++i ){
    const double x = loc_grid_dim.get_sq_ctr_x(i);
    i_on_loc_grid_edge = ((0 == i) || (loc_grid_dim.get_w() == (i+1))) ? 
        true : false;
    for( j = ij_min.second; j <= ij_max.second; ++j ){
        j_on_loc_grid_edge = ((0 == j) || (loc_grid_dim.get_h() == (j+1))) ?
            true : false;
        if(i_on_loc_grid_edge || j_on_loc_grid_edge){
            should_add = true;
            }
        else{
            const double y = loc_grid_dim.get_sq_ctr_y(j);
            const double fwd_dot = (m_fwd.get_x() * x) + (m_fwd.get_y() * y);
            const double fwd_cross = (m_fwd.get_x() * y) - (m_fwd.get_y() * x);
            should_add = ((min_fwd_dot<=fwd_dot) && (fwd_dot<=max_fwd_dot) &&
                (min_fwd_cross<=fwd_cross) && (fwd_cross<=max_fwd_cross) ) ?
                true : false;
            }
        if(should_add){
            index_vec->push_back(np02_uint16_pair(i,j));
            }
        }
    }

AA_DECR_CALL_DEPTH();
}

/* signed distance.  positive => point lies outside rectangle.
  negative => point lies inside rectangle */
double np02_rect::get_distance_from_xy(const np02_xy& xy,
    np02_xy *near_xy) const{
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data_num(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_2);
double d = 0.0;
double dx, dy, dsq;
const double fwd_dot = (m_fwd.get_x() * xy.get_x()) + 
                       (m_fwd.get_y() * xy.get_y());
const double fwd_cross = (m_fwd.get_x() * xy.get_y()) - 
                         (m_fwd.get_y() * xy.get_x());
const double fwd_dot_0 = get_fwd_dot_0();
const double fwd_dot_1 = get_fwd_dot_1();
const double fwd_cross_0 = get_fwd_cross_0();
const double fwd_cross_1 = get_fwd_cross_1();


/*                        :               :
                          :               :      fwd_cross > fwd_cross_1
                          :p01            :p11
           - - - - - - - -*---------------*- - - - - - - -
                          |               |                      
                          |      ctr      |                      
                         h|       +-------|--------->fwd         
                          |               |                      
                          |               |                      
           - - - - - - - -*---------------*- - - - - - - -       
                          :p00     w      :p10
                          :               :      fwd_cross < fwd_cross_0
                          :               :  
      fwd_dot < fwd_dot_0 :               : fwd_dot > fwd_dot_1    */
if((fwd_dot < fwd_dot_0) && (fwd_cross < fwd_cross_0)){
    d = m_p00.get_distance_to(xy);
    if(NULL != near_xy){ *near_xy = m_p00; }
    }
else if((fwd_dot < fwd_dot_0) && (fwd_cross > fwd_cross_1)){
    d = m_p01.get_distance_to(xy);
    if(NULL != near_xy){ *near_xy = m_p01; }
    }
else if((fwd_dot > fwd_dot_1) && (fwd_cross < fwd_cross_0)){
    d = m_p10.get_distance_to(xy);
    if(NULL != near_xy){ *near_xy = m_p10; }
    }
else if((fwd_dot > fwd_dot_1) && (fwd_cross > fwd_cross_1)){
    d = m_p11.get_distance_to(xy);
    if(NULL != near_xy){ *near_xy = m_p11; }
    }
else{
/*                        :               :
                          :d_cross > d_fwd: 1
                       p01:               :p11
           - - - - - - - -*---------------*- - - - - - - -
                          | \    ctr    / |                      
          d_fwd > d_cross |  +----+----+  |--------->fwd         
                          | /           \ |   d_fwd > d_cross                   
           - - - - - - - -*---------------*- - - - - - - -       
                       p00:               :p10
                          :d_cross > d_fwd:
                          :               :                */
    const double d_fwd = ( fwd_dot < m_fwd_dot_ctr ) ?
        (fwd_dot_0 - fwd_dot) : (fwd_dot - fwd_dot_1);
    const double d_cross = ( fwd_cross < m_fwd_cross_ctr ) ?
        (fwd_cross_0 - fwd_cross) : (fwd_cross - fwd_cross_1);
    d = (d_fwd > d_cross) ? d_fwd : d_cross;
    if(NULL != near_xy){
        if( ( fwd_dot_0 <= fwd_dot ) && ( fwd_dot <= fwd_dot_1 ) &&
            ( fwd_cross_0 <= fwd_cross ) && ( fwd_cross <= fwd_cross_1 ) ){
            /* xy is inside rectangle */
            *near_xy = xy;
            }
        else{
            double fd, fc;
            if(d_fwd > d_cross){
                fd = (fwd_dot < m_fwd_dot_ctr ) ? fwd_dot_0 : fwd_dot_1;
                fc = fwd_cross;
                }
            else{
                fd = fwd_dot;
                fc = ( fwd_cross < m_fwd_cross_ctr ) ? fwd_cross_0 : fwd_cross_1;
                }
            near_xy->set_x((fd * m_fwd.get_x()) - (fc * m_fwd.get_y()));
            near_xy->set_y((fd * m_fwd.get_y()) + (fc * m_fwd.get_x()));
            }
        }
    }
AA_DECR_CALL_DEPTH();
return d;
}

/* signed distance.  positive => segment lies outside rectangle.
  negative => segment lies at least partially inside rectangle */
double np02_rect::get_distance_from_line_seg_ab(const np02_xy& xy_a,
    const np02_xy& xy_b, np02_xy *near_xy,
    np02_xy *other_near_xy) const{
AA_INCR_CALL_DEPTH();
/* distance to line segment endpoints */
np02_xy near_xy_a, near_xy_b;
const double d_a = get_distance_from_xy(xy_a, &near_xy_a);
const double d_b = get_distance_from_xy(xy_b, &near_xy_b);
double d = 0.0;
np02_xy oth_near_xy,  nr_xy;
if(d_a < d_b){
    d = d_a;
    nr_xy = near_xy_a;
    oth_near_xy = xy_a;
    }
else{
    d = d_b;
    nr_xy = near_xy_b;
    oth_near_xy = xy_b;
    }

/* check for distance from rectangle corners to line segment. */
const np02_xy dab(xy_b.get_x()-xy_a.get_x(), xy_b.get_y()-xy_a.get_y());
const double dab_len_sq=(dab.get_x()*dab.get_x())+(dab.get_y()*dab.get_y());
if( dab_len_sq > 1e-40){
    const double dab_len=sqrt(dab_len_sq);
    const np02_xy dap00(m_p00.get_x()-xy_a.get_x(),m_p00.get_y()-xy_a.get_y());
    const np02_xy dap01(m_p01.get_x()-xy_a.get_x(),m_p01.get_y()-xy_a.get_y());
    const np02_xy dap10(m_p10.get_x()-xy_a.get_x(),m_p10.get_y()-xy_a.get_y());
    const np02_xy dap11(m_p11.get_x()-xy_a.get_x(),m_p11.get_y()-xy_a.get_y());

    const double dab_dot_dap00=(dab.get_x() * dap00.get_x()) +
                               (dab.get_y() * dap00.get_y());
    if((dab_dot_dap00 >= 0.0) && (dab_dot_dap00 <= dab_len_sq)){
        const double dab_cross_dap00=(dab.get_x() * dap00.get_y()) -
                                     (dab.get_y() * dap00.get_x());
        const double d_p00 = fabs(dab_cross_dap00)/dab_len;
        if(d_p00 < d){
            d = d_p00;
            nr_xy = m_p00;
            const double f_nr_p00 = dab_dot_dap00/dab_len_sq;
            oth_near_xy.set_x(xy_a.get_x() + (dab.get_x() * f_nr_p00));
            oth_near_xy.set_y(xy_a.get_y() + (dab.get_y() * f_nr_p00));
            }
        }

    const double dab_dot_dap01=(dab.get_x() * dap01.get_x()) +
                               (dab.get_y() * dap01.get_y());
    if((dab_dot_dap01 >= 0.0) && (dab_dot_dap01 <= dab_len_sq)){
        const double dab_cross_dap01=(dab.get_x() * dap01.get_y()) -
                                     (dab.get_y() * dap01.get_x());
        const double d_p01 = fabs(dab_cross_dap01)/dab_len;
        if(d_p01 < d){
            d = d_p01;
            nr_xy = m_p01;
            const double f_nr_p01 = dab_dot_dap01/dab_len_sq;
            oth_near_xy.set_x(xy_a.get_x() + (dab.get_x() * f_nr_p01));
            oth_near_xy.set_y(xy_a.get_y() + (dab.get_y() * f_nr_p01));
            }
        }

    const double dab_dot_dap10=(dab.get_x() * dap10.get_x()) +
                               (dab.get_y() * dap10.get_y());
    if((dab_dot_dap10 >= 0.0) && (dab_dot_dap10 <= dab_len_sq)){
        const double dab_cross_dap10=(dab.get_x() * dap10.get_y()) -
                                     (dab.get_y() * dap10.get_x());
        const double d_p10 = fabs(dab_cross_dap10)/dab_len;
        if(d_p10 < d){
            d = d_p10;
            nr_xy = m_p10;
            const double f_nr_p10 = dab_dot_dap10/dab_len_sq;
            oth_near_xy.set_x(xy_a.get_x() + (dab.get_x() * f_nr_p10));
            oth_near_xy.set_y(xy_a.get_y() + (dab.get_y() * f_nr_p10));
            }
        }

    const double dab_dot_dap11=(dab.get_x() * dap11.get_x()) +
                               (dab.get_y() * dap11.get_y());
    if((dab_dot_dap11 >= 0.0) && (dab_dot_dap11 <= dab_len_sq)){
        const double dab_cross_dap11=(dab.get_x() * dap11.get_y()) -
                                     (dab.get_y() * dap11.get_x());
        const double d_p11 = fabs(dab_cross_dap11)/dab_len;
        if(d_p11 < d){
            d = d_p11;
            nr_xy = m_p11;
            const double f_nr_p11 = dab_dot_dap11/dab_len_sq;
            oth_near_xy.set_x(xy_a.get_x() + (dab.get_x() * f_nr_p11));
            oth_near_xy.set_y(xy_a.get_y() + (dab.get_y() * f_nr_p11));
            }
        }
    }

/* check for bounding box overlap */
np02_xy seg_bb_xy_min, seg_bb_xy_max;
if(xy_a.get_x() < xy_b.get_x()){
    seg_bb_xy_min.set_x(xy_a.get_x());
    seg_bb_xy_max.set_x(xy_b.get_x());
    }
else{
    seg_bb_xy_min.set_x(xy_b.get_x());
    seg_bb_xy_max.set_x(xy_a.get_x());
    }
if(xy_a.get_y() < xy_b.get_y()){
    seg_bb_xy_min.set_y(xy_a.get_y());
    seg_bb_xy_max.set_y(xy_b.get_y());
    }
else{
    seg_bb_xy_min.set_y(xy_b.get_y());
    seg_bb_xy_max.set_y(xy_a.get_y());
    }
np02_xy rect_bb_xy_min, rect_bb_xy_max;
get_bb(&rect_bb_xy_min, &rect_bb_xy_max);
bool x_bb_overlap = false;

if( /* rect min X within seg X range */
    ( (seg_bb_xy_min.get_x() <= rect_bb_xy_min.get_x()) && 
      (rect_bb_xy_min.get_x() <= seg_bb_xy_max.get_x()) ) ||

    /* rect max X within seg X range */
    ( (seg_bb_xy_min.get_x() <= rect_bb_xy_max.get_x()) &&
      (rect_bb_xy_max.get_x() <= seg_bb_xy_max.get_x()) ) ||

    /* seg min X within rect X range */
    ( (rect_bb_xy_min.get_x() <= seg_bb_xy_min.get_x()) &&
      (seg_bb_xy_min.get_x() <= rect_bb_xy_max.get_x()) ) ||

    /* seg max X within rect X range */
    ( (rect_bb_xy_min.get_x() <= seg_bb_xy_max.get_x()) &&
      (seg_bb_xy_max.get_x() <= rect_bb_xy_max.get_x()) ) ){
    x_bb_overlap=true;}

bool y_bb_overlap = false;
if( /* rect min Y within seg Y range */
    ( (seg_bb_xy_min.get_y() <= rect_bb_xy_min.get_y()) && 
      (rect_bb_xy_min.get_y() <= seg_bb_xy_max.get_y()) ) ||

    /* rect max Y within seg Y range */
    ( (seg_bb_xy_min.get_y() <= rect_bb_xy_max.get_y()) &&
      (rect_bb_xy_max.get_y() <= seg_bb_xy_max.get_y()) ) ||

    /* seg min Y within rect Y range */
    ( (rect_bb_xy_min.get_y() <= seg_bb_xy_min.get_y()) &&
      (seg_bb_xy_min.get_y() <= rect_bb_xy_max.get_y()) ) ||

    /* seg max Y within rect Y range */
    ( (rect_bb_xy_min.get_y() <= seg_bb_xy_max.get_y()) &&
      (seg_bb_xy_max.get_y() <= rect_bb_xy_max.get_y()) ) ){
    y_bb_overlap=true;}

if(x_bb_overlap && y_bb_overlap){
    /* bounding boxes overlap, so do further distance checks */
    const np02_xy xy_ab_ave((xy_a.get_x()+xy_b.get_x())/2.0,
                                (xy_a.get_y()+xy_b.get_y())/2.0);
    const double dx_ab = xy_b.get_x() - xy_a.get_x();
    const double dy_ab = xy_b.get_y() - xy_a.get_y();
    const double dsq_ab = (dx_ab*dx_ab) + (dy_ab*dy_ab);
    np02_xy fwd_ab;
    if( dsq_ab > 1e-80){
        const double d_ab = sqrt(dsq_ab);
        fwd_ab.set_x(dx_ab/d_ab);
        fwd_ab.set_y(dy_ab/d_ab);
        }
    else{ fwd_ab.set_x(1.0); fwd_ab.set_y(0.0); }
    const double fwd_ab_dot_a = (fwd_ab.get_x() * xy_a.get_x()) +
                                (fwd_ab.get_y() * xy_a.get_y());
    const double fwd_ab_dot_b = (fwd_ab.get_x() * xy_b.get_x()) +
                                (fwd_ab.get_y() * xy_b.get_y());
    const double fwd_ab_cross_ab = (fwd_ab.get_x() * xy_ab_ave.get_y()) -
                                   (fwd_ab.get_y() * xy_ab_ave.get_x());

    /* near points on line seg in perpendicular zone of line seg 
      and within rectangle

                       A
                        \
               +---------\-------+corner
               |          \      |
               | near point*     |
               |            \    |
               |             \   |
               +--------------\--+
                               \
                                \
                                 B
    */
    np02_xy xy_np00, xy_np01, xy_np10, xy_np11; 
    bool np00_in_zone, np01_in_zone, np10_in_zone, np11_in_zone;
    double np00_d, np01_d, np10_d, np11_d;
    for(size_t i = 0; i < 4; ++i){
        const np02_xy *corner_xy;
        np02_xy *xy_np;
        bool *np_in_zone;
        double *np_d;
        switch(i){
            case 0: corner_xy = &m_p00; xy_np=&xy_np00; 
                np_in_zone=&np00_in_zone; np_d = &np00_d; break;
            case 1: corner_xy = &m_p01; xy_np=&xy_np01; 
                np_in_zone=&np01_in_zone; np_d = &np01_d; break;
            case 2: corner_xy = &m_p10; xy_np=&xy_np10; 
                np_in_zone=&np10_in_zone; np_d = &np10_d; break;
            case 3: default: corner_xy = &m_p11; xy_np=&xy_np11; 
                np_in_zone=&np11_in_zone; np_d = &np11_d; break;
            }
        xy_np->set_x(0.0);
        xy_np->set_y(0.0);
        *np_in_zone = false;
        *np_d = 0.0;
        const double fwd_ab_dot_corner = 
            (fwd_ab.get_x() * (corner_xy->get_x())) +
            (fwd_ab.get_y() * (corner_xy->get_y()));
        if((fwd_ab_dot_a <= fwd_ab_dot_corner) &&
           (fwd_ab_dot_corner <= fwd_ab_dot_b)){
            /* corner is in seg perpendicular zone, so near point is as well */
            const np02_xy cross_ab(-fwd_ab.get_y(), fwd_ab.get_x());
            xy_np->set_x((fwd_ab_dot_corner * fwd_ab.get_x()) +
                         (fwd_ab_cross_ab * cross_ab.get_x()));
            xy_np->set_y((fwd_ab_dot_corner * fwd_ab.get_y()) +
                         (fwd_ab_cross_ab * cross_ab.get_y()));

            /* check if xy_near_point is within */
            const double rect_fwd_dot_np = (m_fwd.get_x() * (xy_np->get_x())) +
                                           (m_fwd.get_y() * (xy_np->get_y()));
            const double rect_fwd_cross_np=(m_fwd.get_x() * (xy_np->get_y())) -
                                           (m_fwd.get_y() * (xy_np->get_x()));
            if( (fabs(m_fwd_dot_ctr - rect_fwd_dot_np) < ( m_w/2.0 )) &&
                (fabs(m_fwd_cross_ctr - rect_fwd_cross_np) < ( m_h/2.0 )) ){
                /* near point is within rectangle */
                *np_in_zone = true;
                const double fwd_ab_cross_corner = 
                    (cross_ab.get_x() * (corner_xy->get_x())) +
                    (cross_ab.get_y() * (corner_xy->get_y()));
                *np_d = -fabs(fwd_ab_cross_corner - fwd_ab_cross_ab);

                if(AA_SHOULD_RUN_XDBG(CF01_AA_DEBUG_LEVEL_2)){
                    const double np_d_double_check_dx =
                        xy_np->get_x() - (corner_xy->get_x());
                    const double np_d_double_check_dy =
                        xy_np->get_y() - (corner_xy->get_y());
                    const double np_d_double_check = sqrt(
                        (np_d_double_check_dx*np_d_double_check_dx) +
                        (np_d_double_check_dy*np_d_double_check_dy) );
                    const double np_d_err=fabs(np_d_double_check-fabs(*np_d));
                    const double max_np_d_err = (np_d_double_check > 1.0) ?
                        np_d_double_check * 1e-8 : 1e-8;
                    AA_XDBG_ASSERT(np_d_err <= max_np_d_err,CF01_AA_DEBUG_LEVEL_2 );
                    }
                }
            }
        }

    /*                 A
                        \
               +---------\-------+p11
               |          \      |
               |           *np11 |
               |            \    |
               |         np00*   |
               |              \  |
               |               \ |
            p00+----------------\+
                                 \
                                  \
                                   B
    */  
    if(np00_in_zone && np11_in_zone){
        /* take higher value => signed dist. to closer corner of p00 and p11 */
        if(np00_d < np11_d){
            if(np11_d < d){
                d = np11_d;
                nr_xy = m_p11;
                oth_near_xy = xy_np11;
                }
            }
        else{
            if(np00_d < d){
                d = np00_d;
                nr_xy = m_p00;
                oth_near_xy = xy_np00;
                }
            }
        }
    if(np01_in_zone && np10_in_zone){
        /* take higher value => signed dist. to closer corner of p00 and p11 */
        if(np01_d < np10_d){
            if(np10_d < d){
                d = np10_d;
                nr_xy = m_p10;
                oth_near_xy = xy_np10;
                }
            }
        else{
            if(np01_d < d){
                d = np01_d;
                nr_xy = m_p01;
                oth_near_xy = xy_np01;
                }
            }
        }      
    }

if(NULL != near_xy){ *near_xy = nr_xy; }
if(NULL != other_near_xy){ *other_near_xy = oth_near_xy; }
AA_DECR_CALL_DEPTH();
return d;
}

double np02_rect::get_distance_from_circle(const np02_circle *c,
    np02_xy *near_xy, np02_xy *circle_near_xy) const{
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != c);
AUTO_ASSERT(c->get_radius() >= 0.0);
np02_xy nr_xy;
const double d_ctr = get_distance_from_xy(c->get_ctr(), &nr_xy);
const double d = d_ctr - c->get_radius();
if(NULL != near_xy){*near_xy = nr_xy;}
if(NULL != circle_near_xy){
    if( d > 0.0 ){
        const double& c_ctr_x = (c->get_ctr()).get_x();
        const double& c_ctr_y = (c->get_ctr()).get_y();
        const double& c_radius = c->get_radius();
        const double dx = nr_xy.get_x() - c_ctr_x;
        const double dy = nr_xy.get_y() - c_ctr_y;
        const double ddsq = (dx*dx) + (dy*dy);
        const double dd = sqrt(ddsq);
        AUTO_ASSERT(fabs(dd - fabs(d_ctr)) <= fabs((dd + fabs(d_ctr))*1e-6) );
        np02_xy unit_v(1.0,0.0);
        if( dd < 1e-40){
            if(fabs(dx) > fabs(dy)){
                unit_v.set_x((dx>0.0) ? 1.0:-1.0);
                unit_v.set_y(0.0);
                }
            else{
                unit_v.set_x(0.0);
                unit_v.set_y((dy>0.0) ? 1.0:-1.0);
                }
            }
        else{
            unit_v.set_x(dx/dd);
            unit_v.set_y(dy/dd);
            }
        circle_near_xy->set_x(c_ctr_x + (c_radius * unit_v.get_x()));
        circle_near_xy->set_y(c_ctr_y + (c_radius * unit_v.get_y()));
        }
    else if( d_ctr > 0.0 ){
        *circle_near_xy = nr_xy;
        }
    else{
        *circle_near_xy = c->get_ctr();
        }
    }
AA_DECR_CALL_DEPTH();
return d;
}

double np02_rect::get_distance_from_line_seg(
    const np02_line_seg *n, np02_xy *near_xy,
    np02_xy *line_seg_near_xy)const{
const np02_xy& xy_a = n->get_p_0();
const np02_xy& xy_b = n->get_p_1();
np02_xy nr_xy;
np02_xy line_seg_ctr_xy;
const double d_ctr = get_distance_from_line_seg_ab(xy_a, xy_b, &nr_xy,
    &line_seg_ctr_xy);
double d = d_ctr - ((n->get_width())/2.0);

np02_xy line_seg_near_xy_0;
np02_xy ln_seg_near_xy;
const double hw = ( n->get_width() ) / 2.0;
if( d_ctr > hw ){
    const double dx = line_seg_ctr_xy.get_x() - nr_xy.get_x();
    const double dy = line_seg_ctr_xy.get_y() - nr_xy.get_y();
    const double ddsq = (dx*dx) + (dy*dy);
    const double dd = sqrt(ddsq);
    AUTO_ASSERT(fabs(dd - fabs(d_ctr)) <= fabs((dd + fabs(d_ctr))*1e-6) );
    np02_xy unit_v(1.0,0.0);
    if( dd < 1e-40){
        if(fabs(dx) > fabs(dy)){
            unit_v.set_x((dx>0.0) ? 1.0:-1.0);
            unit_v.set_y(0.0);
            }
        else{
            unit_v.set_x(0.0);
            unit_v.set_y((dy>0.0) ? 1.0:-1.0);
            }
        }
    else{
        unit_v.set_x(dx/dd);
        unit_v.set_y(dy/dd);
        }
    line_seg_near_xy_0.set_x(nr_xy.get_x() + (d * unit_v.get_x()));
    line_seg_near_xy_0.set_y(nr_xy.get_y() + (d * unit_v.get_y()));
    }
else{
    line_seg_near_xy_0 = nr_xy;
    }
n->get_distance_from_xy(line_seg_near_xy_0, &ln_seg_near_xy);


np02_xy nr_xy2, line_seg_near_xy_3;
const double d2 = get_distance_from_xy(ln_seg_near_xy, &nr_xy2);
const double d3 = n->get_distance_from_xy(nr_xy, &line_seg_near_xy_3);
if(d2 < d){
    d = d2;
    nr_xy = nr_xy2;
    }
if(d3 < d){
    d = d3;
    ln_seg_near_xy = line_seg_near_xy_3;
    }

if(NULL != near_xy){ *near_xy = nr_xy; }
if( NULL != line_seg_near_xy){ *line_seg_near_xy = ln_seg_near_xy; }
return d;
}

double np02_rect::get_distance_from_rect(const np02_rect *r,
    np02_xy *near_xy, np02_xy *rect_near_xy) const{
/*

    p01     p11              p11
     *-------*                *
     |\     /|              / | \
     | \   / |            /   |   \
     |  \ /  |       p01*-----+-----*p10
     |   +   |            \   |   /
     |  / \  |              \ | /
     | /   \ |                *
     |/     \|               p00
     *-------*
    p00     p10
*/
np02_xy nr_xy_a, oth_nr_xy_a,
            nr_xy_b, oth_nr_xy_b,
            nr_xy_c, oth_nr_xy_c,
            nr_xy_d, oth_nr_xy_d,
            nr_xy, oth_nr_xy;
const double d_a = r->get_distance_from_line_seg_ab(m_p00, m_p11, &oth_nr_xy_a, &nr_xy_a);
const double d_b = r->get_distance_from_line_seg_ab(m_p10, m_p01, &oth_nr_xy_b, &nr_xy_b);
const double d_c = get_distance_from_line_seg_ab(r->m_p00, r->m_p11, &nr_xy_c, &oth_nr_xy_c);
const double d_d = get_distance_from_line_seg_ab(r->m_p10, r->m_p01, &nr_xy_d, &oth_nr_xy_d);

double d = d_a;
nr_xy = nr_xy_a;
oth_nr_xy = oth_nr_xy_a;
if(d_b < d){ 
    d = d_b;
    nr_xy = nr_xy_b;
    oth_nr_xy = oth_nr_xy_b;
    }
if(d_c < d){ 
    d = d_c;
    nr_xy = nr_xy_c;
    oth_nr_xy = oth_nr_xy_c; 
    }
if(d_d < d){ 
    d = d_d;
    nr_xy = nr_xy_d;
    oth_nr_xy = oth_nr_xy_d;
    }

/*

    p01  H  p11              p11
     *---+---*                *
     |   |   |             H/   \F
     |   |   |            / \   / \
     |   |   |       p01*     +    *p10
    E+---+---+F           \ /   \ /
     |   |   |             E\   /G
     |   |   |                *
     |   |   |               p00
     *---+---*
    p00  G  p10
*/
np02_xy xy_e((m_p00.get_x()+m_p01.get_x())/2.0,(m_p00.get_y()+m_p01.get_y())/2.0);
np02_xy xy_f((m_p10.get_x()+m_p11.get_x())/2.0,(m_p10.get_y()+m_p11.get_y())/2.0);
np02_xy xy_g((m_p00.get_x()+m_p10.get_x())/2.0,(m_p00.get_y()+m_p10.get_y())/2.0);
np02_xy xy_h((m_p01.get_x()+m_p11.get_x())/2.0,(m_p01.get_y()+m_p11.get_y())/2.0);

np02_xy oth_xy_e(((r->m_p00).get_x()+(r->m_p01).get_x())/2.0,((r->m_p00).get_y()+(r->m_p01).get_y())/2.0);
np02_xy oth_xy_f(((r->m_p10).get_x()+(r->m_p11).get_x())/2.0,((r->m_p10).get_y()+(r->m_p11).get_y())/2.0);
np02_xy oth_xy_g(((r->m_p00).get_x()+(r->m_p10).get_x())/2.0,((r->m_p00).get_y()+(r->m_p10).get_y())/2.0);
np02_xy oth_xy_h(((r->m_p01).get_x()+(r->m_p11).get_x())/2.0,((r->m_p01).get_y()+(r->m_p11).get_y())/2.0);

np02_xy nr_xy_5, oth_nr_xy_5,
            nr_xy_6, oth_nr_xy_6,
            nr_xy_7, oth_nr_xy_7,
            nr_xy_8, oth_nr_xy_8;
const double d_5 = r->get_distance_from_line_seg_ab(xy_e, xy_f, &oth_nr_xy_5, &nr_xy_5);
const double d_6 = r->get_distance_from_line_seg_ab(xy_g, xy_h, &oth_nr_xy_6, &nr_xy_6);
const double d_7 = get_distance_from_line_seg_ab(oth_xy_e, oth_xy_f, &nr_xy_7, &oth_nr_xy_7);
const double d_8 = get_distance_from_line_seg_ab(oth_xy_g, oth_xy_h, &nr_xy_8, &oth_nr_xy_8);

if(d_5 < d){ 
    d = d_5;
    nr_xy = nr_xy_5;
    oth_nr_xy = oth_nr_xy_5;
    }
if(d_6 < d){ 
    d = d_6;
    nr_xy = nr_xy_6;
    oth_nr_xy = oth_nr_xy_6;
    }
if(d_7 < d){ 
    d = d_7;
    nr_xy = nr_xy_7;
    oth_nr_xy = oth_nr_xy_7;
    }
if(d_8 < d){ 
    d = d_8;
    nr_xy = nr_xy_8;
    oth_nr_xy = oth_nr_xy_8;
    }


np02_xy  oth_nr_xy_9, nr_xy_10;
const double d_9  = r->get_distance_from_xy(nr_xy, &oth_nr_xy_9);
const double d_10 = get_distance_from_xy(oth_nr_xy, &nr_xy_10);
if(d_9 < d){ 
    d = d_9;
    oth_nr_xy = oth_nr_xy_9;
    }
if(d_10 < d){ 
    d = d_10;
    nr_xy = nr_xy_10;
    }


if(NULL != near_xy){ *near_xy = nr_xy; }
if(NULL != rect_near_xy){ *rect_near_xy = oth_nr_xy; }
return d;
}

double np02_rect::get_distance_from_polygon(const np02_polygon *p,
    np02_xy *near_xy, np02_xy *rect_near_xy) const{
/* TODO: implement */
return 0.0;
}

double np02_rect::get_distance_from_shape(const np02_shape *s,
    np02_xy *near_xy, np02_xy *other_near_xy) const{
double d = 0.0;
AA_INCR_CALL_DEPTH();
if(NULL != s){
    const np02_circle *c = dynamic_cast<const np02_circle *>(s);
    const np02_line_seg *n = dynamic_cast<const np02_line_seg *>(s);
    const np02_rect *r = dynamic_cast<const np02_rect *>(s);
    const np02_polygon *p = dynamic_cast<const np02_polygon *>(s);
    if(NULL != c){ d = get_distance_from_circle(c, near_xy, other_near_xy); }
    else if(NULL!=n){ d=get_distance_from_line_seg(n, near_xy, other_near_xy);}
    else if(NULL!=r){ d=get_distance_from_rect(r, near_xy, other_near_xy); }
    else if(NULL!=p){ d=get_distance_from_polygon(p, near_xy, other_near_xy); }
    else{ AA_ALWAYS_ASSERT(false); }
    }
AA_DECR_CALL_DEPTH();
return d;
}

void np02_rect::translate_no_loc_grid(const np02_xy& dxy){
const np02_xy ctr(m_ctr.get_x() + dxy.get_x(), m_ctr.get_y() + dxy.get_y());
init(ctr, m_w, m_h, m_rot_deg);
}

void np02_rect::rotate_no_loc_grid(const np02_xy& rot_ctr,
    const double& rot_deg){
double cos_rot = 1.0;
double sin_rot = 0.0;
if(0.0 == rot_deg){ cos_rot = 1.0; sin_rot = 0.0; }
else if(90.0 == rot_deg){ cos_rot = 0.0; sin_rot = 1.0; }
else if(-90.0 == rot_deg){ cos_rot = 0.0; sin_rot = -1.0; }
else if(180.0 == rot_deg){ cos_rot = -1.0; sin_rot = 0.0; }
else if(-180.0 == rot_deg){ cos_rot = -1.0; sin_rot = 0.0; }
else{
    const double rot_rad =
        rot_deg * (3.1415926535897932384626433832795 / 180.0);
    cos_rot = cos(rot_rad);
    sin_rot = sin(rot_rad);
    }
const np02_xy rot_arm_initial(m_ctr.get_x() - rot_ctr.get_x(),
    m_ctr.get_y() - rot_ctr.get_y());
const np02_xy rot_arm_final(
    (rot_arm_initial.get_x() * cos_rot) -
    (rot_arm_initial.get_y() * sin_rot),
    (rot_arm_initial.get_x() * sin_rot) +
    (rot_arm_initial.get_y() * cos_rot));
const np02_xy ctr(rot_ctr.get_x() + rot_arm_final.get_x(),
    rot_ctr.get_y() + rot_arm_final.get_y());
double rot_d = m_rot_deg + rot_deg;
if( rot_d > 360.0){ rot_d -= 360.0; }
else if( rot_d < -360.0){ rot_d -= 360.0; }
init(ctr, m_w, m_h, rot_d);
}

int np02_rect::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const{
int err_cnt;
err_cnt += np02_shape::verify_data(err_msg,err_msg_capacity,err_msg_pos);

if(NP02_SHAPE_TYPE_RECT != get_shape_type() ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "rect: this=%x  shape_type=%i\n", this, get_shape_type()); 
    }

const np02_shp_alloc *shp_alloc = get_shp_alloc();
if((NULL != shp_alloc) && 
    (this != shp_alloc->alloc_get_rect_by_idx(get_shape_idx()))){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "circle: this=%x  != (shp_alloc=%x) rect_by_idx(idx=%i)=%x\n",
        this, shp_alloc, get_shape_idx(),
        shp_alloc->alloc_get_rect_by_idx(get_shape_idx())); 
    }

err_cnt += verify_data_num(err_msg,err_msg_capacity,err_msg_pos);
err_cnt += verify_data_rect_loc_grid(err_msg,err_msg_capacity,err_msg_pos);
return 0;
}

int np02_rect::verify_data_num( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const{
int err_cnt = 0;
double dx, dy;
if(m_w < 0.0){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }
if(m_h < 0.0){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }
if(fabs(m_rot_deg) >= 360.0){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }
const double rot_rad = m_rot_deg * (3.1415926535897932384626433832795/180.0);
const double cos_rot = cos(rot_rad);
const double sin_rot = sin(rot_rad);
if(fabs(m_fwd.get_x() - cos_rot) > 1e-8){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }
if(fabs(m_fwd.get_y() - sin_rot) > 1e-8){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }

const double fwd_len_sq = (m_fwd.get_x() * m_fwd.get_x()) +
                          (m_fwd.get_y() * m_fwd.get_y());
if(fabs(fwd_len_sq - 1.0) > 1e-8){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }
const double ctr_len_sq = (m_ctr.get_x() * m_ctr.get_x()) +
                          (m_ctr.get_y() * m_ctr.get_y());
const double p00_len_sq = (m_p00.get_x() * m_p00.get_x()) +
                          (m_p00.get_y() * m_p00.get_y());
const double p01_len_sq = (m_p01.get_x() * m_p01.get_x()) +
                          (m_p01.get_y() * m_p01.get_y());
const double p10_len_sq = (m_p10.get_x() * m_p10.get_x()) +
                          (m_p10.get_y() * m_p10.get_y());
const double p11_len_sq = (m_p11.get_x() * m_p11.get_x()) +
                          (m_p11.get_y() * m_p11.get_y());
const double w_sq = m_w * m_w;
const double h_sq = m_h * m_h;
double len_sq_criteria = ctr_len_sq;
if(p00_len_sq > len_sq_criteria){len_sq_criteria = p00_len_sq;}
if(p01_len_sq > len_sq_criteria){len_sq_criteria = p01_len_sq;}
if(p10_len_sq > len_sq_criteria){len_sq_criteria = p10_len_sq;}
if(p11_len_sq > len_sq_criteria){len_sq_criteria = p11_len_sq;}
if(w_sq > len_sq_criteria){len_sq_criteria = w_sq;}
if(h_sq > len_sq_criteria){len_sq_criteria = h_sq;}
const double max_d_sq_err =
    (len_sq_criteria > 1.0) ? len_sq_criteria * 1e-8 : 1e-8;
const double max_d_err =
    (len_sq_criteria > 1.0) ? sqrt(len_sq_criteria) * 1e-8 : 1e-8;

const double fwd_dot_ctr =(m_fwd.get_x() * m_ctr.get_x()) +
                          (m_fwd.get_y() * m_ctr.get_y());
const double fwd_cross_ctr =(m_fwd.get_x() * m_ctr.get_y()) -
                            (m_fwd.get_y() * m_ctr.get_x());
if(fabs(fwd_dot_ctr - m_fwd_dot_ctr) > max_d_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }
if(fabs(fwd_cross_ctr - m_fwd_cross_ctr) > max_d_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }

dx = m_p10.get_x() - m_p00.get_x();
dy = m_p10.get_y() - m_p00.get_y();
const double w0_sq = (dx*dx) + (dy*dy);
if( fabs(w0_sq - w_sq) > max_d_sq_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s  dx=%g  dy=%g  w0_sq=%g  w_sq=%g  "
        "max_d_sq_err=%g\n", __FILE__, __LINE__, __FUNCTION__,
        dx, dy, w0_sq, w_sq, max_d_sq_err );
    }
if(fabs((m_fwd.get_x()*m_w) - dx) > max_d_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }
if(fabs((m_fwd.get_y()*m_w) - dy) > max_d_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }

dx = m_p11.get_x() - m_p01.get_x();
dy = m_p11.get_y() - m_p01.get_y();
const double w1_sq = (dx*dx) + (dy*dy);
if( fabs(w1_sq - w_sq) > max_d_sq_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s  dx=%g  dy=%g  w1_sq=%g  w_sq=%g  "
        "max_d_sq_err=%g\n", __FILE__, __LINE__, __FUNCTION__,
        dx, dy, w1_sq, w_sq, max_d_sq_err );
    }
if(fabs((m_fwd.get_x()*m_w) - dx) > max_d_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }
if(fabs((m_fwd.get_y()*m_w) - dy) > max_d_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }

dx = m_p01.get_x() - m_p00.get_x();
dy = m_p01.get_y() - m_p00.get_y();
const double h0_sq = (dx*dx) + (dy*dy);
if( fabs(h0_sq - h_sq) > max_d_sq_err){
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s  dx=%g  dy=%g  h0_sq=%g  h_sq=%g  "
        "max_d_sq_err=%g\n", __FILE__, __LINE__, __FUNCTION__,
        dx, dy, h0_sq, h_sq, max_d_sq_err );
    ++err_cnt;
    }
if(fabs((-m_fwd.get_y()*m_h) - dx) > max_d_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }
if(fabs((m_fwd.get_x()*m_h) - dy) > max_d_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }

dx = m_p11.get_x() - m_p10.get_x();
dy = m_p11.get_y() - m_p10.get_y();
const double h1_sq = (dx*dx) + (dy*dy);
if( fabs(h1_sq - h_sq) > max_d_sq_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s  dx=%g  dy=%g  h1_sq=%g  h_sq=%g  "
        "max_d_sq_err=%g\n", __FILE__, __LINE__, __FUNCTION__,
        dx, dy, h1_sq, h_sq, max_d_sq_err );
    }
if(fabs((-m_fwd.get_y()*m_h) - dx) > max_d_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }
if(fabs((m_fwd.get_x()*m_h) - dy) > max_d_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }

const double fwd_dot_0 = get_fwd_dot_0();
const double fwd_dot_00 = (m_fwd.get_x() * m_p00.get_x()) +
                          (m_fwd.get_y() * m_p00.get_y());
const double fwd_dot_01 = (m_fwd.get_x() * m_p01.get_x()) +
                          (m_fwd.get_y() * m_p01.get_y());
if(fabs(fwd_dot_0 - fwd_dot_00) > max_d_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }
if(fabs(fwd_dot_0 - fwd_dot_01) > max_d_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }

const double fwd_dot_1 = get_fwd_dot_1();
const double fwd_dot_10 = (m_fwd.get_x() * m_p10.get_x()) +
                          (m_fwd.get_y() * m_p10.get_y());
const double fwd_dot_11 = (m_fwd.get_x() * m_p11.get_x()) +
                          (m_fwd.get_y() * m_p11.get_y());
if(fabs(fwd_dot_1 - fwd_dot_10) > max_d_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }
if(fabs(fwd_dot_1 - fwd_dot_11) > max_d_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }

const double fwd_cross_0 = get_fwd_cross_0();
const double fwd_cross_00 = (m_fwd.get_x() * m_p00.get_y()) -
                            (m_fwd.get_y() * m_p00.get_x());
const double fwd_cross_10 = (m_fwd.get_x() * m_p10.get_y()) -
                            (m_fwd.get_y() * m_p10.get_x());
if(fabs(fwd_cross_0 - fwd_cross_00) > max_d_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }
if(fabs(fwd_cross_0 - fwd_cross_10) > max_d_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }

const double fwd_cross_1 = get_fwd_cross_1();
const double fwd_cross_01 = (m_fwd.get_x() * m_p01.get_y()) -
                            (m_fwd.get_y() * m_p01.get_x());
const double fwd_cross_11 = (m_fwd.get_x() * m_p11.get_y()) -
                            (m_fwd.get_y() * m_p11.get_x());
if(fabs(fwd_cross_1 - fwd_cross_01) > max_d_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }
if(fabs(fwd_cross_1 - fwd_cross_11) > max_d_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }
return err_cnt;
}

int np02_rect::verify_data_rect_loc_grid( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;
const np02_loc_grid_node *lg_node = get_head_loc_grid_node();
const np02_loc_grid *loc_grid = (NULL == lg_node) ?
    NULL : lg_node->get_loc_grid();
const double extra_search_d = (NULL == loc_grid) ?
    0.0 : loc_grid->get_extra_search_d();
static const np02_loc_grid_dim default_loc_grid_dim;
const np02_loc_grid_dim& loc_grid_dim = (NULL == loc_grid) ?
    default_loc_grid_dim : loc_grid->get_loc_grid_dim();
const size_t max_loc_grid_node_count = 
    static_cast<size_t>(loc_grid_dim.get_w()) *
    static_cast<size_t>(loc_grid_dim.get_h());
/* TODO: check that each locator grid square is within
    (extra_search_d + sqrt(2)* sq_sz) from from rectangle */

/* TODO: check that, if head locator grid is not NULL, that
all needed grid squares are on the list */

return err_cnt;
}

std::ostream& np02_rect::ostream_output(std::ostream& os) const{
os << "<rect>\n";
os << std::hex;
os << "<this>" << this << "</this>\n";
os << std::dec;
np02_shape::ostream_output(os);
os << "<ctr><x>" << m_ctr.get_x() << "</x>"
     << "<y>" << m_ctr.get_y() << "</y></ctr>\n";
os << "<w>" << m_w << "</w>\n";
os << "<h>" << m_h << "</h>\n";
os << "<rot_deg>" << m_rot_deg << "</rot_deg>\n";
os << "<fwd><x>" << m_fwd.get_x() << "</x>"
     << "<y>" << m_fwd.get_y() << "</y></fwd>\n";
os << "<fwd_dot_ctr>" << m_fwd_dot_ctr << "</fwd_dot_ctr>\n";
os << "<fwd_cross_ctr>" << m_fwd_cross_ctr << "</fwd_cross_ctr>\n";
os << "<p00><x>" << m_p00.get_x() << "</x>"
     << "<y>" << m_p00.get_y() << "</y></p00>\n";
os << "<p01><x>" << m_p01.get_x() << "</x>"
     << "<y>" << m_p01.get_y() << "</y></p01>\n";
os << "<p10><x>" << m_p10.get_x() << "</x>"
     << "<y>" << m_p10.get_y() << "</y></p10>\n";
os << "<p11><x>" << m_p11.get_x() << "</x>"
     << "<y>" << m_p11.get_y() << "</y></p11>\n";
os << "</rect>\n";

return os;
}

void np02_rect::write_bmp_file(const np02_xy& xy_min,
    const double& pixel_num, const np02_bmp_color& color,
    np02_bmp_file *bmp_file) const{
np02_shape::write_bmp_file(xy_min, pixel_num, color, bmp_file);
}

void np02_rect::write_dxf_file(const std::string& layer,
    const uint8_t& color, np02_dxf_file *dxf_file) const{
AA_INCR_CALL_DEPTH();
AUTO_ASSERT(NULL != dxf_file);
dxf_file->draw_line(layer,m_p00.get_x(),m_p00.get_y(),m_p01.get_x(),m_p01.get_y(),color);
dxf_file->draw_line(layer,m_p01.get_x(),m_p01.get_y(),m_p11.get_x(),m_p11.get_y(),color);
dxf_file->draw_line(layer,m_p11.get_x(),m_p11.get_y(),m_p10.get_x(),m_p10.get_y(),color);
dxf_file->draw_line(layer,m_p10.get_x(),m_p10.get_y(),m_p00.get_x(),m_p00.get_y(),color);
AA_DECR_CALL_DEPTH();
}


np02_line_seg::np02_line_seg():m_p_0(0.0,0.0), m_p_1(0.0,0.0),
    m_width(0.0),m_fwd(1.0,0.0),m_fwd_dot_0(0.0), m_fwd_dot_1(0.0),
    m_fwd_cross_01(0.0){
set_shape_type(NP02_SHAPE_TYPE_LINE_SEG);
}

np02_line_seg::~np02_line_seg(){}

void np02_line_seg::init(const np02_xy& p_0,const np02_xy& p_1,
    const double& width){
AA_INCR_CALL_DEPTH();
AUTO_ASSERT(width >= 0.0);
m_p_0 = p_0; 
m_p_1 = p_1;
m_width = width;
const double dx = p_1.get_x() - p_0.get_x();
const double dy = p_1.get_y() - p_0.get_y();
const double dsq = (dx*dx) + (dy*dy);
if(dsq < 1e-80){
    if(fabs(dx) > fabs(dy)){
        m_fwd.set_x((dx>0.0)?1.0:-1.0);
        m_fwd.set_y(0.0);
        }
    else{
        m_fwd.set_x(0.0);
        m_fwd.set_y((dy>0.0)?1.0:-1.0);
        }
    }
else{
    const double d = sqrt(dsq);
    m_fwd.set_x(dx/d);
    m_fwd.set_y(dy/d);
    }
const np02_xy ctr01((p_0.get_x() + p_1.get_x())/2.0,
                        (p_0.get_y() + p_1.get_y())/2.0);
m_fwd_dot_0 = (m_fwd.get_x() * p_0.get_x()) + (m_fwd.get_y() * p_0.get_y());
m_fwd_dot_1 = (m_fwd.get_x() * p_1.get_x()) + (m_fwd.get_y() * p_1.get_y());
m_fwd_cross_01=(m_fwd.get_x()*ctr01.get_y()) - (m_fwd.get_y()*ctr01.get_x());

AA_XDBG_ASSERT(0 == verify_data_num(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_1);
AA_DECR_CALL_DEPTH();
}

void np02_line_seg::get_bb(np02_xy *xy_min, np02_xy *xy_max) const{
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != xy_min);
AA_ALWAYS_ASSERT(NULL != xy_max);
AUTO_ASSERT(m_width >= 0.0);
const double hw = m_width/2.0;
if( m_p_0.get_x() < m_p_1.get_x() ){
    xy_min->set_x( m_p_0.get_x() - hw );
    xy_max->set_x( m_p_1.get_x() + hw ); 
    }
else{
    xy_min->set_x( m_p_1.get_x() - hw );
    xy_max->set_x( m_p_0.get_x() + hw ); 
    }
if( m_p_0.get_y() < m_p_1.get_y() ){
    xy_min->set_y( m_p_0.get_y() - hw );
    xy_max->set_y( m_p_1.get_y() + hw ); 
    }
else{
    xy_min->set_y( m_p_1.get_y() - hw );
    xy_max->set_y( m_p_0.get_y() + hw ); 
    }
AA_DECR_CALL_DEPTH();
}

void np02_line_seg::get_loc_grid_indices_for_init(
    const np02_loc_grid_dim& loc_grid_dim,
    const double& extra_search_d, np02_uint16_pair_vec *index_vec)
    const{
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != index_vec);
AUTO_ASSERT(m_width >= 0.0);
AUTO_ASSERT(extra_search_d >= 0.0);
AA_XDBG_ASSERT(0 == loc_grid_dim.verify_data(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_2);

uint16_t i,j;
np02_xy xy_min, xy_max;
get_bb(&xy_min, &xy_max);
xy_min.set_x(xy_min.get_x()-extra_search_d);
xy_min.set_y(xy_min.get_y()-extra_search_d);
xy_max.set_x(xy_max.get_x()+extra_search_d);
xy_max.set_y(xy_max.get_y()+extra_search_d);

np02_uint16_pair ij_min, ij_max;
loc_grid_dim.get_bb_indices(xy_min, xy_max, &ij_min, &ij_max);

const double d_threshold = (m_width/2.0) + extra_search_d +
    (0.75 * loc_grid_dim.get_sq_size());
const double dsq_threshold = d_threshold * d_threshold;
AUTO_ASSERT(d_threshold > 0.0);

const double min_fwd_dot = m_fwd_dot_0 - d_threshold;
const double max_fwd_dot = m_fwd_dot_1 + d_threshold;
const double min_fwd_cross = m_fwd_cross_01 - d_threshold;
const double max_fwd_cross = m_fwd_cross_01 + d_threshold;

bool i_on_loc_grid_edge, j_on_loc_grid_edge, should_add;
for( i = ij_min.first; i <= ij_max.first; ++i ){
    const double x = loc_grid_dim.get_sq_ctr_x(i);
    i_on_loc_grid_edge = ((0 == i) || (loc_grid_dim.get_w() == (i+1))) ? 
        true : false;
    for( j = ij_min.second; j <= ij_max.second; ++j ){
        j_on_loc_grid_edge = ((0 == j) || (loc_grid_dim.get_h() == (j+1))) ?
            true : false;
        if(i_on_loc_grid_edge || j_on_loc_grid_edge){
            should_add = true;
            }
        else{
            const double y = loc_grid_dim.get_sq_ctr_y(j);
            const double fwd_dot = (m_fwd.get_x() * x) + (m_fwd.get_y() * y);
            const double fwd_cross = (m_fwd.get_x() * y) - (m_fwd.get_y() * x);
            should_add = ( (min_fwd_dot<=fwd_dot) && (fwd_dot<=max_fwd_dot) &&
                (min_fwd_cross<=fwd_cross) && (fwd_cross<=max_fwd_cross) ) ?
                true : false;
            }
        if(should_add){
            index_vec->push_back(np02_uint16_pair(i,j));
            }
        }
    }
AA_DECR_CALL_DEPTH();
}

double np02_line_seg::get_distance_from_xy(const np02_xy& xy,
    np02_xy *near_xy) const{
AA_INCR_CALL_DEPTH();
AUTO_ASSERT(m_width >= 0.0);
AA_XDBG_ASSERT(0 == verify_data_num(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_2);

double dx, dy, dsq;
double d_ctr = 0.0;
const double hw = m_width / 2.0;
const double fwd_dot = ( m_fwd.get_x() * xy.get_x() ) +
                       ( m_fwd.get_y() * xy.get_y() );
if(fwd_dot < m_fwd_dot_0){
    /* p_0 is closer to xy than any point on line segment */
    dx = xy.get_x() - m_p_0.get_x();
    dy = xy.get_y() - m_p_0.get_y();
    dsq = (dx*dx) + (dy*dy);
    d_ctr = sqrt(dsq);
    if(NULL != near_xy){
        if( d_ctr > hw ){
            near_xy->set_x(m_p_0.get_x() + ((hw * dx) / d_ctr)); 
            near_xy->set_y(m_p_0.get_y() + ((hw * dy) / d_ctr)); 
            }
        else{
            *near_xy = xy;
            }
        }
    }
else if(fwd_dot > m_fwd_dot_1){
    /* p_1 is closer to xy than any point on line segment */
    dx = xy.get_x() - m_p_1.get_x();
    dy = xy.get_y() - m_p_1.get_y();
    dsq = (dx*dx) + (dy*dy);
    d_ctr = sqrt(dsq);
    if(NULL != near_xy){
        if( d_ctr > hw ){ 
            near_xy->set_x(m_p_1.get_x() + ((hw * dx) / d_ctr)); 
            near_xy->set_y(m_p_1.get_y() + ((hw * dy) / d_ctr)); 
            }
        else{
            *near_xy = xy;
            }
        }
    }
else{
    /* near point lies along line segment between p_0 and p_1  */
    const double fwd_cross = ( m_fwd.get_x() * xy.get_y() ) -
                             ( m_fwd.get_y() * xy.get_x() );
    d_ctr = fabs(fwd_cross - m_fwd_cross_01);
    if(NULL != near_xy){
        if( d_ctr > hw ){
            const double cross = m_fwd_cross_01 + ((fwd_cross > m_fwd_cross_01) ? hw : -hw );
            near_xy->set_x( (fwd_dot * m_fwd.get_x()) - (cross * m_fwd.get_y()) );
            near_xy->set_y( (fwd_dot * m_fwd.get_y()) + (cross * m_fwd.get_x()) );
            }
        else{
            *near_xy = xy;
            }
        }
    }

const double d = d_ctr - hw;
AA_DECR_CALL_DEPTH();
return d;
}

double np02_line_seg::get_distance_from_line_seg_ab(const np02_xy& xy_a,
    const np02_xy& xy_b, np02_xy *near_xy,
    np02_xy *other_near_xy) const{
/* TODO: implement */
return 0.0;
}

double np02_line_seg::get_distance_from_circle(
    const np02_circle *c, np02_xy *near_xy,
    np02_xy *circle_near_xy) const{
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != c);
AUTO_ASSERT(c->get_radius() >= 0.0);
AUTO_ASSERT(m_width >= 0.0);
AA_XDBG_ASSERT(0 == verify_data_num(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_2);
const double d = c->get_distance_from_line_seg(this, circle_near_xy, near_xy);
AA_DECR_CALL_DEPTH();
return d;
}

double np02_line_seg::get_distance_from_line_seg(
    const np02_line_seg *n, np02_xy *near_xy,
    np02_xy *line_seg_near_xy)const{
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != n);
AUTO_ASSERT(n->get_width() >= 0.0);
AUTO_ASSERT(m_width >= 0.0);
AA_XDBG_ASSERT(0 == verify_data_num(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_2);
AA_XDBG_ASSERT(0 == n->verify_data_num(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_2);

const double mhw = m_width / 2.0;
const double nhw = (n->get_width()) / 2.0;
np02_xy xya, xyb, xyc, xyd;
const double d_m0 = get_distance_from_xy(n->get_p_0(),&xya) - nhw;
const double d_m1 = get_distance_from_xy(n->get_p_1(),&xyb) - nhw;
const double d_n0 = n->get_distance_from_xy(get_p_0(),&xyc) - mhw;
const double d_n1 = n->get_distance_from_xy(get_p_1(),&xyd) - mhw;
double d = 0.0;
np02_xy nr_xy(0.0,0.0);
np02_xy oth_nr_xy(0.0,0.0);

if( (d_m0 < d_m1) && (d_m0 < d_n0) && (d_m0 < d_n1)){
    d = d_m0;
    nr_xy = xya;
    n->get_distance_from_xy(xya, &oth_nr_xy);
    }
else if( (d_m1 < d_n0) && (d_m1 < d_n1)){
    d = d_m1;
    nr_xy = xyb;
    n->get_distance_from_xy(xyb, &oth_nr_xy);
    } 
else if(d_n0 < d_n1){
    d = d_n0;
    oth_nr_xy = xyc;
    get_distance_from_xy(xyc, &nr_xy);
    } 
else{
    d = d_n1;
    oth_nr_xy = xyd;
    get_distance_from_xy(xyd, &nr_xy);
    }

/* check for intersection: if quadrilateral (M0,N1,M1,N0) is convex 

         N0......M1
        .  *   /   .
        .    /*     .
       .   /     *   .
       . /          * .
      M0...............N1
*/
const double& xm0 = get_p_0().get_x();
const double& ym0 = get_p_0().get_y();
const double& xm1 = get_p_1().get_x();
const double& ym1 = get_p_1().get_y();
const double& xn0 = (n->get_p_0()).get_x();
const double& yn0 = (n->get_p_0()).get_y();
const double& xn1 = (n->get_p_1()).get_x();
const double& yn1 = (n->get_p_1()).get_y();
const double dxm0n1 = xn1 - xm0;
const double dym0n1 = yn1 - ym0;
const double dxn1m1 = xm1 - xn1;
const double dyn1m1 = ym1 - yn1;
const double dxm1n0 = xn0 - xm1;
const double dym1n0 = yn0 - ym1;
const double dxn0m0 = xm0 - xn0;
const double dyn0m0 = ym0 - yn0;
const double cross_m0n1m1 = (dxm0n1*dyn1m1)-(dym0n1*dxn1m1);
const double cross_n1m1n0 = (dxn1m1*dym1n0)-(dyn1m1*dxm1n0);
const double cross_m1n0m0 = (dxm1n0*dyn0m0)-(dym1n0*dxn0m0);
const double cross_n0m0n1 = (dxn0m0*dym0n1)-(dyn0m0*dxm0n1);
const bool intersect = ( ( ( cross_m0n1m1 <= 0.0) &&
                           ( cross_n1m1n0 <= 0.0) &&
                           ( cross_m1n0m0 <= 0.0) &&
                           ( cross_n0m0n1 <= 0.0) ) ||
                         ( ( cross_m0n1m1 >= 0.0) &&
                           ( cross_n1m1n0 >= 0.0) &&
                           ( cross_m1n0m0 >= 0.0) &&
                           ( cross_n0m0n1 >= 0.0) ) );
if(intersect){
    /*    
         N0......M1
        .  * P /   .
        .    /*     .
       .   /     *len_n
       . /len_m     * .
      M0...............N1
    */
    const double len_m = m_fwd_dot_1-m_fwd_dot_0;
    const double len_n = (n->get_fwd_dot_1())-(n->get_fwd_dot_0());
    if((len_m > 1e-40) && (len_n > 1e-40)){
        /* Distance away is negative because one line segment
        must be moved that distance just to get to zero. 
        Distance is measured perpendicular from segment to
        endpoint of other segment
                 N0......M1               N0......M1      
                .  * dm1/  .             .  \dn0*   .     
                .   / *     .            .    *  \   .    
               .  /      *   .          .   *      \  .   
               ./dm0        * .         . *       dn1\ .  
              M0...............N1      M0...............N1
        */
        const double dn1 = -fabs(cross_m0n1m1/len_m);
        const double dm1 = -fabs(cross_n1m1n0/len_n);
        const double dn0 = -fabs(cross_m1n0m0/len_m);
        const double dm0 = -fabs(cross_n0m0n1/len_n);

        /* chose shortest abs distance to any */
        double d_ctr = dn1;
        if(dm1 > d_ctr){d_ctr = dm1;}
        if(dn0 > d_ctr){d_ctr = dn0;}
        if(dn1 > d_ctr){d_ctr = dn1;}
        const double dd = d_ctr - (mhw + nhw);

        if(dd < d){
            d = dd; 

            /*
            find intersection point (px, py)
            cross products provide 2 equations, 2 unknowns
            [-m_fwd.y   m_fwd.x ] [px]  =  m_fwd_cross_01;
            [-n_fwd.y   n_fwd.x ] [py]  =  n_fwd_cross_01;

            inv([a b]   = (1/det)[ d -b ], det=ad-bc
                [c d])           [-c  a ]
            */
            const double& m_fwd_x = m_fwd.get_x();
            const double& m_fwd_y = m_fwd.get_y();
            const double& n_fwd_x = (n->get_fwd()).get_x();
            const double& n_fwd_y = (n->get_fwd()).get_y();
            const double& n_fwd_cross_01 = n->get_fwd_cross_01();
            const double det = (-m_fwd_y * n_fwd_x) + ( m_fwd_x * n_fwd_y);
            if( fabs(det) > 1e-20){
                nr_xy.set_x((( n_fwd_x * m_fwd_cross_01 ) - 
                             ( m_fwd_x * n_fwd_cross_01 ))/det );
                nr_xy.set_y((( n_fwd_y * m_fwd_cross_01 ) - 
                             ( m_fwd_y * n_fwd_cross_01 ))/det );
                }
            else{
                /* take the average of all endpoints */
                nr_xy.set_x((xm0 + xm1 + xn0 + xn1)/4.0 );
                nr_xy.set_y((ym0 + ym1 + yn0 + yn1)/4.0 );
                }
            oth_nr_xy = nr_xy;
            }
        }
    }


if( NULL != near_xy ){ *near_xy = nr_xy; } 
if( NULL != line_seg_near_xy ){ *line_seg_near_xy = oth_nr_xy; } 

if(AA_SHOULD_RUN_XDBG(CF01_AA_DEBUG_LEVEL_1)){
    double d2_err_estimate = 0.0;
    np02_xy near_xy2, line_seg_near_xy2;
    double d2 = get_distance_from_line_seg_double_check(n, &d2_err_estimate,
        &near_xy2, &line_seg_near_xy2);
    const double near_xy_err = nr_xy.get_distance_to(near_xy2);
    const double line_seg_near_xy_err =
        oth_nr_xy.get_distance_to(line_seg_near_xy2);

    if( d > 0.0 ){
        /*
        if(fabs(d - d2) >= (2 * d2_err_estimate)){
            std::cout << "fabs(d - d2) >= (2 * d2_err_estimate)\n";
            std::cout << "d = " << d << "\n";
            std::cout << "d2 = " << d2 << "\n";
            std::cout << "d2_err_estimate = " << d2_err_estimate << "\n";
            std::cout << "nr_xy = " << nr_xy.get_x() << "," << nr_xy.get_y() << "\n";
            std::cout << "near_xy2 = " << near_xy2.get_x() << "," << near_xy2.get_y() << "\n";
            std::cout << "oth_nr_xy = " << oth_nr_xy.get_x() << "," << oth_nr_xy.get_y() << "\n";
            std::cout << "line_seg_near_xy2 = " << line_seg_near_xy2.get_x()
                << "," << line_seg_near_xy2.get_y() << "\n";
            }
        */

        AA_XDBG_ASSERT(fabs(d - d2) < (2 * d2_err_estimate),CF01_AA_DEBUG_LEVEL_1 );
        AA_XDBG_ASSERT(near_xy_err < (2 * d2_err_estimate),CF01_AA_DEBUG_LEVEL_1 );
        AA_XDBG_ASSERT(line_seg_near_xy_err < (2 * d2_err_estimate),CF01_AA_DEBUG_LEVEL_1 );
        }
    }

AA_DECR_CALL_DEPTH();
return d;
}


double np02_line_seg::get_distance_from_line_seg_double_check(
    const np02_line_seg *n, double *err_estimate,
    np02_xy *near_xy, np02_xy *line_seg_near_xy)const{
AA_INCR_CALL_DEPTH();
double d = 0.0;
np02_xy nr_xy(0.0, 0.0);
np02_xy oth_nr_xy(0.0, 0.0);
if(NULL != n){
    double err_est = 0.0;
    double d_near_m = std::numeric_limits<double>::max();
    double f_near_m = 0.5;
    double f_low_m = 0.0;
    double f_high_m = 1.0;
    np02_xy line_seg_nr_xy_m(0.0, 0.0);
    double d_near_n = std::numeric_limits<double>::max();
    double f_near_n = 0.5;
    double f_low_n = 0.0;
    double f_high_n = 1.0;
    np02_xy line_seg_nr_xy_n(0.0, 0.0);

    for(int i = 0; i < 5; ++i){
        for(int pt_m_or_n = 0; pt_m_or_n < 2; ++pt_m_or_n){
            const np02_line_seg *line_seg;
            const np02_line_seg *p_line_seg;
            double *d_near;
            double *f_near;
            double *f_low;
            double *f_high;
            np02_xy *line_seg_nr_xy;
            switch(pt_m_or_n){
                case 0:
                    line_seg = n;
                    p_line_seg = this;
                    d_near = &d_near_m;
                    f_near = &f_near_m;
                    f_low = &f_low_m;
                    f_high = &f_high_m;
                    line_seg_nr_xy = &line_seg_nr_xy_m;
                    break;
                case 1: /* fall through */
                default:
                    line_seg = this;
                    p_line_seg = n;
                    d_near = &d_near_n;
                    f_near = &f_near_n;
                    f_low = &f_low_n;
                    f_high = &f_high_n;
                    line_seg_nr_xy = &line_seg_nr_xy_n;
                    break;
                }
            const np02_xy& p0 = p_line_seg->get_p_0();
            const np02_xy& p1 = p_line_seg->get_p_1();
            AUTO_ASSERT(*f_low >= 0.0);
            AUTO_ASSERT(*f_low < *f_high);
            AUTO_ASSERT(*f_high <= 1.0);
            const double f_step = (*f_high - *f_low)/4.0;
            const double f_end = *f_low + (5.0 * f_step);
            for( double f = *f_low; f < f_end; f += f_step){
                const double g = 1.0 - f;
                const np02_xy p((g*p0.get_x()) + (f*p1.get_x()),
                                    (g*p0.get_y()) + (f*p1.get_y()));
                np02_xy ln_seg_nr_xy(0.0, 0.0);
                const double d_p =
                    line_seg->get_distance_from_xy(p, &ln_seg_nr_xy) -
                    ((p_line_seg->get_width())/2.0);
                if(d_p < *d_near){
                    *d_near = d_p;
                    *f_near = f;
                    *line_seg_nr_xy = ln_seg_nr_xy;
                    }
                }
            *f_low = *f_near - (f_step / 2.0);
            if(*f_low < 0.0){ *f_low = 0.0; }
            *f_high = *f_near + (f_step / 2.0);
            if(*f_high > 1.0){ *f_high = 1.0; }
            AUTO_ASSERT(*f_low < *f_high);
            }
        }


    if( d_near_m < d_near_n ){
        d = d_near_m;
        err_est = (f_high_m - f_low_m) * (m_fwd_dot_1 - m_fwd_dot_0);
        oth_nr_xy = line_seg_nr_xy_m;
        this->get_distance_from_xy(oth_nr_xy, &nr_xy);
        }
    else{
        d = d_near_n;
        err_est=(f_high_n-f_low_n)*((n->get_fwd_dot_1())-(n->get_fwd_dot_0()));
        nr_xy = line_seg_nr_xy_n;
        n->get_distance_from_xy(nr_xy, &oth_nr_xy);
        }
    if((d < 0.0) && (err_est < fabs(d))){
        err_est = fabs(d);
        }
    if(NULL != err_estimate){ *err_estimate = err_est; }
    }
if(NULL != near_xy){ *near_xy = nr_xy; }
if(NULL != line_seg_near_xy){ *line_seg_near_xy = oth_nr_xy; }
AA_DECR_CALL_DEPTH();
return d;
}

double np02_line_seg::get_distance_from_rect(const np02_rect *r,
    np02_xy *near_xy, np02_xy *rect_near_xy) const{
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != r);
return r->get_distance_from_line_seg(this, rect_near_xy, near_xy);
AA_DECR_CALL_DEPTH();
}

double np02_line_seg::get_distance_from_polygon(
    const np02_polygon *p,  np02_xy *near_xy,
    np02_xy *rect_near_xy) const{
/* TODO: implement */
return 0.0;
}

double np02_line_seg::get_distance_from_shape(const np02_shape *s,
    np02_xy *near_xy, np02_xy *other_near_xy) const{
double d = 0.0;
AA_INCR_CALL_DEPTH();
if(NULL != s){
    const np02_circle *c = dynamic_cast<const np02_circle *>(s);
    const np02_line_seg *n = dynamic_cast<const np02_line_seg *>(s);
    const np02_rect *r = dynamic_cast<const np02_rect *>(s);
    const np02_polygon *p = dynamic_cast<const np02_polygon *>(s);
    if(NULL != c){ d = get_distance_from_circle(c, near_xy, other_near_xy); }
    else if(NULL != n){ d = get_distance_from_line_seg(n, near_xy, other_near_xy); }
    else if(NULL != r){ d = get_distance_from_rect(r, near_xy, other_near_xy); }
    else if(NULL != p){ d = get_distance_from_polygon(p, near_xy, other_near_xy); }
    else{ AA_ALWAYS_ASSERT(false); }
    }
AA_DECR_CALL_DEPTH();
return d;
}

void np02_line_seg::translate_no_loc_grid(const np02_xy& dxy){
const np02_xy p0(m_p_0.get_x() + dxy.get_x(), m_p_0.get_y() + dxy.get_y());
const np02_xy p1(m_p_1.get_x() + dxy.get_x(), m_p_1.get_y() + dxy.get_y());
init(p0, p1, m_width);
}

void np02_line_seg::rotate_no_loc_grid(const np02_xy& rot_ctr,
const double& rot_deg){
double cos_rot = 1.0;
double sin_rot = 0.0;
if(0.0 == rot_deg){ cos_rot = 1.0; sin_rot = 0.0; }
else if(90.0 == rot_deg){ cos_rot = 0.0; sin_rot = 1.0; }
else if(-90.0 == rot_deg){ cos_rot = 0.0; sin_rot = -1.0; }
else if(180.0 == rot_deg){ cos_rot = -1.0; sin_rot = 0.0; }
else if(-180.0 == rot_deg){ cos_rot = -1.0; sin_rot = 0.0; }
else{
    const double rot_rad =
        rot_deg * (3.1415926535897932384626433832795 / 180.0);
    cos_rot = cos(rot_rad);
    sin_rot = sin(rot_rad);
    }
const np02_xy rot_arm_initial_0(
    m_p_0.get_x() - rot_ctr.get_x(),
    m_p_0.get_y() - rot_ctr.get_y());
const np02_xy rot_arm_initial_1(
    m_p_1.get_x() - rot_ctr.get_x(),
    m_p_1.get_y() - rot_ctr.get_y());
const np02_xy rot_arm_final_0(
    (rot_arm_initial_0.get_x() * cos_rot) -
    (rot_arm_initial_0.get_y() * sin_rot),
    (rot_arm_initial_0.get_x() * sin_rot) +
    (rot_arm_initial_0.get_y() * cos_rot));
const np02_xy rot_arm_final_1(
    (rot_arm_initial_1.get_x() * cos_rot) -
    (rot_arm_initial_1.get_y() * sin_rot),
    (rot_arm_initial_1.get_x() * sin_rot) +
    (rot_arm_initial_1.get_y() * cos_rot));
const np02_xy p0(
    m_p_0.get_x() + rot_arm_final_0.get_x(),
    m_p_0.get_y() + rot_arm_final_0.get_y());
const np02_xy p1(
    m_p_1.get_x() + rot_arm_final_1.get_x(),
    m_p_1.get_y() + rot_arm_final_1.get_y());
init(p0, p1, m_width);
}

int np02_line_seg::verify_data( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;
err_cnt += np02_shape::verify_data(err_msg,err_msg_capacity,err_msg_pos);

if(NP02_SHAPE_TYPE_LINE_SEG != get_shape_type() ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "line_seg: this=%x  shape_type=%i\n", this, get_shape_type()); 
    }

const np02_shp_alloc *shp_alloc = get_shp_alloc();
if((NULL != shp_alloc) && 
    (this != shp_alloc->alloc_get_line_seg_by_idx(get_shape_idx()))){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "line_seg: this=%x  != (shp_alloc=%x) line_seg_by_idx(idx=%i)=%x\n",
        this, shp_alloc, get_shape_idx(),
        shp_alloc->alloc_get_line_seg_by_idx(get_shape_idx())); 
    }

err_cnt+=verify_data_num(err_msg,err_msg_capacity,err_msg_pos);
/* TODO: check locator grid nodes */

return err_cnt;
}

int np02_line_seg::verify_data_num( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;

const double dx01 = m_p_1.get_x() - m_p_0.get_x();
const double dy01 = m_p_1.get_y() - m_p_0.get_y();
const double d01 = sqrt((dx01*dx01) + (dy01*dy01));
const double max_d01_err = (d01 > 1.0) ? d01 * 1e-8 : 1e-8;
const double d01_from_dot = m_fwd_dot_1 - m_fwd_dot_0;
const double d01_err = fabs(d01 - d01_from_dot);
if(d01_err > max_d01_err){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "np02_line_seg (%f,%f)-(%f,%f) w=%f  d01=%f != d01_from_dot=%f",
        m_p_0.get_x(), m_p_0.get_y(), m_p_1.get_x(), m_p_1.get_y(),
        m_width, d01, d01_from_dot );
    }
if(d01_from_dot < 0.0){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }

if(d01 > 1e-8){
    const np02_xy fwd(dx01/d01, dy01/d01);
    const double fwd_dot_0 = (fwd.get_x() * m_p_0.get_x()) +
                             (fwd.get_y() * m_p_0.get_y());
    const double fwd_dot_1 = (fwd.get_x() * m_p_1.get_x()) +
                             (fwd.get_y() * m_p_1.get_y());
    const double fwd_cross_0 = (fwd.get_x() * m_p_0.get_y()) -
                               (fwd.get_y() * m_p_0.get_x());
    const double fwd_cross_1 = (fwd.get_x() * m_p_1.get_y()) -
                               (fwd.get_y() * m_p_1.get_x());
    double dot_cross_err_threshold = (1e-9) * (
        fabs(fwd_dot_0) + fabs(fwd_dot_1) + 
        fabs(fwd_cross_0) + fabs(fwd_cross_1) +
        fabs(m_fwd_dot_0) + fabs(m_fwd_dot_1) + 
        fabs(m_fwd_cross_01) );
    if(dot_cross_err_threshold < 1e-8){ dot_cross_err_threshold = 1e-8; }
    if(fabs(m_fwd.get_x()-fwd.get_x()) > 1e-8){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
        }
    if(fabs(m_fwd.get_y()-fwd.get_y()) > 1e-8){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
        }
    if(fabs(m_fwd_dot_0-fwd_dot_0) > dot_cross_err_threshold){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
        }
    if(fabs(m_fwd_dot_1-fwd_dot_1) > dot_cross_err_threshold){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
        }
    if(fabs(m_fwd_cross_01-fwd_cross_0) > dot_cross_err_threshold ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
        }
    if(fabs(m_fwd_cross_01-fwd_cross_1) > dot_cross_err_threshold){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
        }
    }

const np02_xy p_0_from_dot_cross(
    (m_fwd_dot_0 * m_fwd.get_x()) - (m_fwd_cross_01 * m_fwd.get_y()),
    (m_fwd_dot_0 * m_fwd.get_y()) + (m_fwd_cross_01 * m_fwd.get_x()) );
const np02_xy p_1_from_dot_cross(
    (m_fwd_dot_1 * m_fwd.get_x()) - (m_fwd_cross_01 * m_fwd.get_y()),
    (m_fwd_dot_1 * m_fwd.get_y()) + (m_fwd_cross_01 * m_fwd.get_x()) );
if( fabs(m_p_0.get_x() - p_0_from_dot_cross.get_x()) > max_d01_err ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }
if( fabs(m_p_0.get_y() - p_0_from_dot_cross.get_y()) > max_d01_err ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }
if( fabs(m_p_1.get_x() - p_1_from_dot_cross.get_x()) > max_d01_err ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }
if( fabs(m_p_1.get_y() - p_1_from_dot_cross.get_y()) > max_d01_err ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }

return err_cnt;
}

std::ostream& np02_line_seg::ostream_output(std::ostream& os) const{
os << "<line_seg>\n";
os << std::hex;
os << "<this>" << this << "</this>\n";
os << std::dec;
np02_shape::ostream_output(os);
os << "<p_0><x>" << m_p_0.get_x() << "</x>"
    << "<y>" << m_p_0.get_y() << "</y></p_0>\n";
os << "<p_1><x>" << m_p_1.get_x() << "</x>"
    << "<y>" << m_p_1.get_y() << "</y></p_1>\n";
os << "<width>" << m_width << "</width>\n";
os << "<fwd><x>" << m_fwd.get_x() << "</x>"
    << "<y>" << m_fwd.get_y() << "</y></fwd>\n";
os << "<fwd_dot_0>" << m_fwd_dot_0 << "</fwd_dot_0>\n";
os << "<fwd_dot_1>" << m_fwd_dot_1 << "</fwd_dot_1>\n";
os << "<fwd_cross_01>" << m_fwd_cross_01 << "</fwd_cross_01>\n";
os << "</line_seg>\n";
return os;
}

void np02_line_seg::write_bmp_file(const np02_xy& xy_min,
    const double& pixel_num, const np02_bmp_color& color,
    np02_bmp_file *bmp_file) const{
np02_shape::write_bmp_file(xy_min, pixel_num, color, bmp_file);
}

void np02_line_seg::write_dxf_file(const std::string& layer,
    const uint8_t& color, np02_dxf_file *dxf_file) const{
AA_INCR_CALL_DEPTH();
AUTO_ASSERT(NULL != dxf_file);
/* centerline */
dxf_file->draw_line(layer,m_p_0.get_x(),m_p_0.get_y(),
                          m_p_1.get_x(),m_p_1.get_y(),color);

/* boundary */
if( m_width > 0.0 ){
    const double w2 = m_width / 2.0;
    dxf_file->draw_circle(layer,m_p_0.get_x(),m_p_0.get_y(),w2,color);
    dxf_file->draw_circle(layer,m_p_1.get_x(),m_p_1.get_y(),w2,color);
    const np02_xy left_w2( -m_fwd.get_y() * w2, m_fwd.get_x() * w2);
    const np02_xy p0_left( m_p_0.get_x() + left_w2.get_x(),
                               m_p_0.get_y() + left_w2.get_y() );
    const np02_xy p0_right(m_p_0.get_x() - left_w2.get_x(),
                               m_p_0.get_y() - left_w2.get_y() );
    const np02_xy p1_left( m_p_1.get_x() + left_w2.get_x(),
                               m_p_1.get_y() + left_w2.get_y() );
    const np02_xy p1_right(m_p_1.get_x() - left_w2.get_x(),
                               m_p_1.get_y() - left_w2.get_y() );
    dxf_file->draw_line(layer,p0_left.get_x(),p0_left.get_y(),
                              p1_left.get_x(),p1_left.get_y(),color);
    dxf_file->draw_line(layer,p0_right.get_x(),p0_right.get_y(),
                              p1_right.get_x(),p1_right.get_y(),color);
    }
AA_DECR_CALL_DEPTH();
}


np02_polygon::np02_polygon(): m_vertices(){
set_shape_type(NP02_SHAPE_TYPE_POLYGON);
}

np02_polygon::~np02_polygon(){}

void np02_polygon::init(const np02_xy_vec& vertices){
m_vertices = vertices; 
np02_xy_itr unique_end = std::unique(m_vertices.begin(), m_vertices.end());
m_vertices.erase(unique_end, m_vertices.end());
}

void np02_polygon::np02_polygon::get_bb(np02_xy *xy_min,
    np02_xy *xy_max) const{
/* TODO: implement */
}

void np02_polygon::get_loc_grid_indices_for_init(
        const np02_loc_grid_dim& loc_grid_dim,
        const double& extra_search_d, np02_uint16_pair_vec *index_vec)
        const{
/* TODO: implement */
}

double np02_polygon::get_distance_from_xy(const np02_xy& xy,
    np02_xy *near_xy) const{
/* TODO: implement */
return 0.0;
}

double np02_polygon::get_distance_from_line_seg_ab(const np02_xy& xy_a,
    const np02_xy& xy_b, np02_xy *near_xy,
    np02_xy *other_near_xy) const{
/* TODO: implement */
return 0.0;
}

double np02_polygon::get_distance_from_circle(
    const np02_circle *c,
    np02_xy *near_xy, np02_xy *circle_near_xy) const{
/* TODO: implement */
return 0.0;
}

double np02_polygon::get_distance_from_line_seg(
    const np02_line_seg *n, np02_xy *near_xy,
    np02_xy *line_seg_near_xy)const{
/* TODO: implement */
return 0.0;
}

double np02_polygon::get_distance_from_rect(const np02_rect *r,
    np02_xy *near_xy, np02_xy *rect_near_xy) const{
/* TODO: implement */
return 0.0;
}

double np02_polygon::get_distance_from_polygon(
    const np02_polygon *p, np02_xy *near_xy,
    np02_xy *rect_near_xy) const{
/* TODO: implement */
return 0.0;
}

double np02_polygon::get_distance_from_shape(const np02_shape *s,
    np02_xy *near_xy, np02_xy *other_near_xy) const{
/* TODO: implement */
return 0.0;
}

void np02_polygon::translate_no_loc_grid(const np02_xy& dxy){ /* TODO: implement */ }
void np02_polygon::rotate_no_loc_grid(const np02_xy& rot_ctr,  const double& rot_deg){ /* TODO: implement */ }

int np02_polygon::verify_data( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;
err_cnt += np02_shape::verify_data(err_msg, err_msg_capacity, err_msg_pos);

if(NP02_SHAPE_TYPE_POLYGON != get_shape_type() ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "polygon: this=%x  shape_type=%i\n", this, get_shape_type()); 
    }

const np02_shp_alloc *shp_alloc = get_shp_alloc();
if((NULL != shp_alloc) && 
    (this != shp_alloc->alloc_get_polygon_by_idx(get_shape_idx()))){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "polygon: this=%x  != (shp_alloc=%x) polygon_by_idx(idx=%i)=%x\n",
        this, shp_alloc, get_shape_idx(),
        shp_alloc->alloc_get_polygon_by_idx(get_shape_idx())); 
    }



/* TODO: implement */
return err_cnt;
}

std::ostream& np02_polygon::ostream_output(std::ostream& os) const{
return os;
}

void np02_polygon::write_bmp_file(const np02_xy& xy_min,
    const double& pixel_num, const np02_bmp_color& color,
    np02_bmp_file *bmp_file) const{
np02_shape::write_bmp_file(xy_min, pixel_num, color, bmp_file);
}

void np02_polygon::write_dxf_file(const std::string& layer,
    const uint8_t& color, np02_dxf_file *dxf_file) const{

}


void np02_loc_grid_dim::init(
    const np02_loc_grid_dim_init_params *init_params){
AA_INCR_CALL_DEPTH();
if(NULL == init_params){
    AA_ALWAYS_ASSERT(false);
    reset();
    }
else if((init_params->point_count == 0) ||
    (init_params->loc_grid_density <= 0.0) ||
    (init_params->max_loc_grid_sq_count == 0) ||
    ((init_params->bb_min_xy).get_x() > (init_params->bb_max_xy).get_x()) ||
    ((init_params->bb_min_xy).get_y() > (init_params->bb_max_xy).get_y())){
    AA_ALWAYS_ASSERT(init_params->point_count > 0);
    AA_ALWAYS_ASSERT(init_params->loc_grid_density > 0.0);
    AA_ALWAYS_ASSERT(init_params->max_loc_grid_sq_count > 0);
    AA_ALWAYS_ASSERT((init_params->bb_min_xy).get_x() <= 
                     (init_params->bb_max_xy).get_x());
    AA_ALWAYS_ASSERT((init_params->bb_min_xy).get_y() <= 
                     (init_params->bb_max_xy).get_y());
    reset();
    }
else{
    /* bounding box */
    const double& x_min = (init_params->bb_min_xy).get_x();
    const double& x_max = (init_params->bb_max_xy).get_x();
    const double& y_min = (init_params->bb_min_xy).get_y();
    const double& y_max = (init_params->bb_max_xy).get_y();
    const double bb_dx = x_max-x_min;
    const double bb_dy = y_max-y_min;
    const double bb_area = bb_dx*bb_dy;
    AA_ALWAYS_ASSERT(bb_dx >= 0.0);
    AA_ALWAYS_ASSERT(bb_dy >= 0.0);
    AA_ALWAYS_ASSERT(bb_area >= 0.0);

    /* estimate number of grid squares */
    double grid_sq_count_float =
       static_cast<double>(init_params->point_count)/
       (init_params->loc_grid_density);
    const double max_grid_sq_count_float =
        static_cast<double>(init_params->max_loc_grid_sq_count);
    if(grid_sq_count_float > max_grid_sq_count_float){
        grid_sq_count_float = max_grid_sq_count_float; }
    if( grid_sq_count_float < 1.0 ){
        grid_sq_count_float = 1.0; }

    /* estimate grid square size */
    double sq_sz_approx = 1.0;
    if( bb_area > 0.0 ){
        const double sq_area_approx = bb_area/grid_sq_count_float;
        AA_ALWAYS_ASSERT(sq_area_approx > 0.0);
        sq_sz_approx = sqrt(sq_area_approx);
        }
    else if( bb_dx > 0.0 ){
        sq_sz_approx = bb_dx / grid_sq_count_float;
        }
    else if( bb_dy > 0.0 ){
        sq_sz_approx = bb_dy / grid_sq_count_float;
        }
    else{
        sq_sz_approx = 1.0;
        }

    /* round up to an integer/(2^k) fraction */
    int n, f33;
    double f, ff;
    f = frexp( sq_sz_approx, &n );
    AA_ALWAYS_ASSERT(f >= 0.5);
    f33 = 1 + static_cast<int>(32.0 * f);
    ff = static_cast<double>(f33)/32.0;
    AA_ALWAYS_ASSERT(ff <= 1.0);
    AA_ALWAYS_ASSERT(ff >= 0.5);
    set_sq_size(ldexp(ff, n));
    AA_ALWAYS_ASSERT(get_sq_size() >= sq_sz_approx);
    AA_ALWAYS_ASSERT(get_sq_size() > 0.0);

    /* loc grid width, height*/
    static const double max_wh_float = static_cast<double>(
        std::numeric_limits<uint16_t>::max() - 1);
    double w_float = 1.0 + (bb_dx / get_sq_size());
    double h_float = 1.0 + (bb_dy / get_sq_size());
    if( w_float > max_wh_float ){ w_float = max_wh_float; }
    if( h_float > max_wh_float ){ h_float = max_wh_float; }
    set_w( static_cast<uint16_t>(w_float));
    set_h( static_cast<uint16_t>(h_float));

    /* loc grid lower left corner*/
    const double bb_ctr_x = (x_min+x_max)/2.0;
    const double bb_ctr_y = (y_min+y_max)/2.0;
    const double gw = static_cast<double>(get_w())*get_sq_size();
    const double gh = static_cast<double>(get_h())*get_sq_size();
    AUTO_ASSERT(gw >= bb_dx);
    AUTO_ASSERT(gh >= bb_dy);
    set_x_min(bb_ctr_x-(gw/2.0));
    set_y_min(bb_ctr_y-(gh/2.0));
    AUTO_ASSERT(get_x_min() <= x_min);
    AUTO_ASSERT(get_y_min() <= y_min);
    AUTO_ASSERT((get_x_min()+gw) >= x_max);
    AUTO_ASSERT((get_y_min()+gh) >= y_max);
    }
AUTO_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_DECR_CALL_DEPTH();
}

void np02_loc_grid_dim::get_bb_indices(
    const np02_xy& xy_min, const np02_xy& xy_max,
    np02_uint16_pair *ij_min, np02_uint16_pair *ij_max) const{
if( NULL != ij_min ) {
    ij_min->first = get_i(xy_min.get_x());
    ij_min->second = get_j(xy_min.get_y());}
if( NULL != ij_max ) {
    ij_max->first = get_i(xy_max.get_x());
    ij_max->second = get_j(xy_max.get_y());}
}

int np02_loc_grid_dim::verify_data( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;
if(((m_w==0) && (m_h!=0)) || ((m_w!=0) && (m_h==0))){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "w=%i != h=%i\n", m_w, m_h); }
if((m_sq_size < 0.0) || ((m_w!=0) && (m_h!=0) && (m_sq_size <= 0.0))){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "w=%i  h=%i  sq_sz=%f\n", m_w, m_h, m_sq_size); }
return err_cnt;
}

std::ostream& np02_loc_grid_dim::ostream_output(std::ostream& os) const{
os << "<loc_grid_dim>\n";
os << "<w>" << m_w << "</w><h>" << m_h << "</h>\n";
os << "<x_min>" << m_x_min << "</x_min><y_min>" << m_y_min << "</y_min>\n";
os << "<sq_size>" << m_sq_size << "</sq_size>\n";
os << "</loc_grid_dim>\n";
return os;
}

void np02_loc_grid_dim::write_bmp_file(const np02_xy& xy_min,
    const double& pixel_num, const np02_bmp_color& color,
    np02_bmp_file *bmp_file) const{
if(NULL != bmp_file){
    uint16_t w_i, h_j;
    int32_t i0, j0, i1, j1;
    for(w_i = 0; w_i <= m_w; ++w_i ){
        i0 = static_cast<int32_t>( pixel_num *
            (m_x_min + (static_cast<double>(w_i)*m_sq_size) - xy_min.get_x()));
        j0 = static_cast<int32_t>( pixel_num * (m_y_min - xy_min.get_x()) );
        j1 = static_cast<int32_t>( pixel_num *
            (m_y_min + (static_cast<double>(m_h)*m_sq_size) - xy_min.get_x()));
        bmp_file->draw_line(i0, j0, i0, j1, color);
        }
    for(h_j = 0; h_j <= m_h; ++h_j ){
        i0 = static_cast<int32_t>( pixel_num * (m_x_min - xy_min.get_x()));
        i1 = static_cast<int32_t>( pixel_num *
            (m_x_min + (static_cast<double>(m_w)*m_sq_size) - xy_min.get_x()));
        j0 = static_cast<int32_t>( pixel_num *
            (m_y_min + (static_cast<double>(h_j)*m_sq_size) - xy_min.get_x()));
        bmp_file->draw_line(i0, j0, i1, j0, color);
        }
    }
}

np02_shp_alloc *np02_loc_grid_node::get_shp_alloc() const{
return (NULL == m_loc_grid) ? NULL : m_loc_grid->get_shp_alloc();
}

lyr_idx_type np02_loc_grid_node::get_lyr_idx() const{
return (NULL == m_loc_grid) ?
    NP02_LYR_INVALID_IDX : m_loc_grid->get_lyr_idx();
}

int np02_loc_grid_node::verify_data( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;
size_t loc_grid_node_s_count=0;
size_t max_loc_grid_node_s_count=0;
size_t loc_grid_node_count=0;
size_t max_loc_grid_node_count = 0xFFFFFF;
const np02_loc_grid_node *loc_grid_node=NULL;
size_t found_count=0;
np02_loc_grid_dim loc_grid_dim;

const np02_shp_alloc *loc_grid_shp_alloc =
    (NULL == m_loc_grid) ? NULL : m_loc_grid->get_shp_alloc();

if(NULL != loc_grid_shp_alloc){
    if(this != loc_grid_shp_alloc->alloc_get_loc_grid_node_by_idx(m_alloc_idx)){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "loc_grid_node:%x != (shp_alloc=%x)->get node by (idx=%i) = %x\n",
            this, loc_grid_shp_alloc, m_alloc_idx,
            loc_grid_shp_alloc->alloc_get_loc_grid_node_by_idx(m_alloc_idx) );
        }
    }

if( NULL != m_loc_grid ) {
    loc_grid_dim = m_loc_grid->get_loc_grid_dim();
    max_loc_grid_node_s_count = static_cast<size_t>(loc_grid_dim.get_w()) 
        * static_cast<size_t>(loc_grid_dim.get_h());
    if(NULL == loc_grid_shp_alloc){
        max_loc_grid_node_count = 0xFFFFFF;
        }
    else{
        max_loc_grid_node_count = (loc_grid_shp_alloc->alloc_get_circle_count()) +
            (loc_grid_shp_alloc->alloc_get_line_seg_count()) +
            (loc_grid_shp_alloc->alloc_get_rect_count()) +
            (loc_grid_shp_alloc->alloc_get_polygon_count());
        }
    }

if(NULL != m_owner){
    const shp_owner_idx_type shp_owner_idx = m_owner->get_shp_owner_index();
    const np02_shp_alloc *shape_shp_alloc = (NULL == m_owner) ?
        NULL : m_owner->get_shp_alloc();
    if(shape_shp_alloc != loc_grid_shp_alloc){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "loc_grid_node:%x (owner=%x/type=%i)->shp_alloc=%x"
             " != (loc_grid=%x)->shp_alloc=%x\n",
            this, m_owner, m_owner->get_shape_type(),
            shape_shp_alloc, m_loc_grid, loc_grid_shp_alloc);
        }

    /* m_owner->loc_grid_head->s_next->s_next ... this*/
    loc_grid_node = m_owner->get_head_loc_grid_node();
    loc_grid_node_s_count = 0;
    found_count = 0;
    while((NULL != loc_grid_node) &&
        (loc_grid_node_s_count <= max_loc_grid_node_s_count)){
        if(this == loc_grid_node){
            ++found_count; }
        ++loc_grid_node_s_count;
        loc_grid_node = loc_grid_node->get_s_next();
        }
    if(1 != found_count){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "loc_grid_node:%x m_owner->loc_grid_head->"
            "s_next found_count=%i\n", this, found_count);
        }
    }

if(NULL != m_loc_grid && NULL != m_owner){
    if(m_i >= loc_grid_dim.get_w()){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "loc_grid_node:%x  i=%i >= w=%i\n", this, m_i,
            loc_grid_dim.get_w());
        }
    if(m_j >= loc_grid_dim.get_h()){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "loc_grid_node:%x  j=%i >= h=%i\n", this, m_j,
            loc_grid_dim.get_h());
        }

    /* m_loc_grid->get_loc_grid_head(m_i,m_j)->next->next-> ... this */
    loc_grid_node = m_loc_grid->get_loc_grid_head_node(m_i, m_j);
    loc_grid_node_count = 0;
    found_count = 0;
    while((NULL != loc_grid_node) &&
        (loc_grid_node_count <= max_loc_grid_node_count)){
        if(this == loc_grid_node){
            ++found_count; }
        ++loc_grid_node_count;
        loc_grid_node = loc_grid_node->get_next();
        }
    if(1 != found_count){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "loc_grid_node:%x loc_grid->loc_grid_head(i=%i,j=%i)->next"
            "->next->... found_count=%i\n", this, m_i, m_j, found_count);
        }
    }

if(NULL != m_loc_grid && NULL == m_owner){
    /* should be found on shp_alloc->m_loc_grid_node_free_chain */
    }
    
if(NULL != m_prev){
    if(this != m_prev->m_next){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "loc_grid_node: this=%x != (prev=%x)->next=%x\n",
            this, m_prev, m_prev->m_next);
        }
    if((m_prev->m_owner == m_owner) && (NULL != m_owner)){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "loc_grid_node=%x: (prev=%x)->owner == owner=%x\n",
            this, m_prev, m_owner);
        }
    if(m_prev->m_loc_grid != m_loc_grid ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "loc_grid_node=%x: (prev=%x)->loc_grid=%x != loc_grid=%x\n",
            this, m_prev, m_prev->m_loc_grid, m_loc_grid);
        }
    }
             
if(NULL != m_next){
    if(this != m_next->m_prev){
        if( NULL == m_loc_grid ){
            /* ok m_next = free chain next */
            }
        else{
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "loc_grid_node: this=%x != (next=%x)->prev=%x\n",
                this, m_next, m_next->m_prev);
            }
        }
    if((m_next->m_owner == m_owner) && (NULL != m_owner) ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "loc_grid_node=%x: (next=%x)->owner == owner=%x\n",
            this, m_next, m_owner);
        }
    if(m_next->m_loc_grid != m_loc_grid ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "loc_grid_node=%x: (next=%x)->loc_grid=%x != loc_grid=%x\n",
            this, m_next, m_next->m_loc_grid, m_loc_grid);
        }
    }
    
if(NULL != m_s_prev){
    if(this != m_s_prev->m_s_next){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "loc_grid_node: this=%x != (s_prev=%x)->s_next=%x\n",
            this, m_s_prev, m_s_prev->m_s_next);
        }
    if(m_s_prev->m_owner != m_owner ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "loc_grid_node=%x: (s_prev=%x)->owner=%x != owner=%x\n",
            this, m_s_prev, m_s_prev->m_owner, m_owner);
        }
    if(m_s_prev->m_loc_grid != m_loc_grid ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "loc_grid_node=%x: (s_prev=%x)->loc_grid=%x != loc_grid=%x\n",
            this, m_s_prev, m_s_prev->m_loc_grid, m_loc_grid);
        }
    }

if(NULL != m_s_next){
    if(this != m_s_next->m_s_prev){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "loc_grid_node: this=%x != (s_next=%x)->s_prev=%x\n",
            this, m_s_next, m_s_next->m_s_prev);
        }
    if(m_s_next->m_owner != m_owner ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "loc_grid_node=%x: (s_next=%x)->owner=%x != owner=%x\n",
            this, m_s_next, m_s_next->m_owner, m_owner);
        }
    if(m_s_next->m_loc_grid != m_loc_grid ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "loc_grid_node=%x: (s_next=%x)->loc_grid=%x != loc_grid=%x\n",
            this, m_s_next, m_s_next->m_loc_grid, m_loc_grid);
        }
    }
return err_cnt;
}

std::ostream& np02_loc_grid_node::ostream_output(std::ostream& os) const{
os << "<loc_grid_node=" << std::hex << this << ">\n";
os << "<owner>" << m_owner << "</owner>\n";
os << "<loc_grid>" << m_loc_grid << "</loc_grid>\n";
os << std::dec << "<i>" << m_i << "</i><j>" << m_j << "</j>\n";
os << std::hex << "<prev>" << m_prev << "</prev>"
               << "<next>" << m_next << "</next>\n";
os << "<s_prev>" << m_s_prev << "</s_prev>"
   << "<s_next>" << m_s_next << "</s_next>\n";
os << std::dec << "</loc_grid_node>\n";
return os;
}

void np02_loc_grid_node::write_bmp_file(const np02_xy& xy_min,
    const double& pixel_num, const np02_bmp_color& color,
    np02_bmp_file *bmp_file) const{

}


np02_loc_grid::np02_loc_grid():
    m_shp_alloc(NULL), m_alloc_idx(0), m_lyr_idx(NP02_LYR_INVALID_IDX),
    m_loc_grid_dim(), m_extra_search_d(0.0), m_loc_grid_vec(),
    m_idx_shape_vec(), m_idx_pair_vec()
{}

np02_loc_grid::~np02_loc_grid(){
AA_INCR_CALL_DEPTH();
clear_loc_grid();
AA_DECR_CALL_DEPTH();
}

void np02_loc_grid::init_loc_grid( const np02_loc_grid_dim& d ){
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_ALWAYS_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()));
size_t loc_grid_sz;
if(!m_loc_grid_vec.empty()){
    clear_loc_grid();}
AUTO_ASSERT(m_loc_grid_vec.empty());
m_loc_grid_dim = d;
loc_grid_sz = static_cast<size_t>(d.get_w()) * static_cast<size_t>(d.get_h());
m_loc_grid_vec.resize(loc_grid_sz, NULL);
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}

void np02_loc_grid::insert_shape_in_loc_grid(np02_shape *shape){
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != shape);
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AUTO_ASSERT(0 == shape->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AUTO_ASSERT(!m_loc_grid_vec.empty());
AUTO_ASSERT(m_loc_grid_dim.get_w() > 0);
AUTO_ASSERT(m_loc_grid_dim.get_h() > 0);
AUTO_ASSERT(m_loc_grid_dim.get_sq_size() > 0.0);
if(NULL != shape){
    /* remove from existing grid */
    np02_loc_grid_node *existing_head_node =
        shape->get_head_loc_grid_node();
    if(NULL != existing_head_node){
        np02_loc_grid *existing_loc_grid =
            existing_head_node->get_loc_grid();
        AA_ALWAYS_ASSERT(NULL != existing_loc_grid);
        if(NULL != existing_loc_grid){
            existing_loc_grid->remove_shape_from_loc_grid(shape);
            }
        }
    AUTO_ASSERT(NULL == shape->get_head_loc_grid_node());

    /* get list of (i,j) indices for shape */
    AUTO_ASSERT(m_idx_pair_vec.empty());
    shape->get_loc_grid_indices_for_init( m_loc_grid_dim, m_extra_search_d,
        &m_idx_pair_vec);

    /* Iterate (i,j) list.  Add node in each grid square */
    np02_uint16_pair_vec_citr ij_itr = m_idx_pair_vec.begin();
    for(; ij_itr != m_idx_pair_vec.end(); ++ij_itr){
        const uint16_t& i = ij_itr->first;
        const uint16_t& j = ij_itr->second;
        const size_t loc_grid_vec_idx = static_cast<size_t>(j) + 
          (static_cast<size_t>(i)*static_cast<size_t>(m_loc_grid_dim.get_h()));
        np02_loc_grid_node *lg_next = m_loc_grid_vec.at(loc_grid_vec_idx);
        np02_loc_grid_node *lg_node = alloc_loc_grid_node();
        lg_node->set_owner(shape);
        AUTO_ASSERT(lg_node->get_loc_grid() == this);
        lg_node->set_i(i);
        lg_node->set_j(j);
        AUTO_ASSERT(lg_node->get_prev() == NULL);
        lg_node->set_next(lg_next);
        if(NULL != lg_next){
            lg_next->set_prev(lg_node); }
        m_loc_grid_vec[loc_grid_vec_idx] = lg_node;
        AUTO_ASSERT(lg_node->get_s_prev() == NULL);
        np02_loc_grid_node *lg_s_next =
            shape->get_head_loc_grid_node();
        if(NULL != lg_s_next){
            lg_s_next->set_s_prev(lg_node); }
        lg_node->set_s_next(lg_s_next);
        shape->set_head_loc_grid_node(lg_node);
        }
    m_idx_pair_vec.clear();
    AUTO_ASSERT(NULL != shape->get_head_loc_grid_node());
    AUTO_ASSERT(0 == shape->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
        AA_ERR_BUF_POS_PTR() ));
    }
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}

void np02_loc_grid::remove_shape_from_loc_grid(np02_shape *shape){
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_ALWAYS_ASSERT(NULL != shape);
AUTO_ASSERT(0 == shape->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AUTO_ASSERT(this == shape->get_loc_grid());
np02_loc_grid_node *lg_node = shape->get_head_loc_grid_node();
while(NULL != lg_node){
    np02_loc_grid_node *lg_prev = lg_node->get_prev();
    np02_loc_grid_node *lg_next = lg_node->get_next();
    np02_loc_grid_node *lg_s_next = lg_node->get_s_next();
    if(NULL == lg_prev){
        const size_t lg_idx = static_cast<size_t>(lg_node->get_j()) +
            (static_cast<size_t>(lg_node->get_i())* 
             static_cast<size_t>(m_loc_grid_dim.get_h()));

        AA_ALWAYS_ASSERT(lg_idx < m_loc_grid_vec.size());
        if(lg_idx < m_loc_grid_vec.size()){
            AUTO_ASSERT(m_loc_grid_vec.at(lg_idx) == lg_node);
            m_loc_grid_vec[lg_idx] = lg_next;
            }
        }
    else{
        lg_prev->set_next(lg_next);
        }
    if(NULL != lg_next){
        lg_next->set_prev(lg_prev);
        }
    free_loc_grid_node(lg_node);
    lg_node=lg_s_next;
    }
shape->set_head_loc_grid_node(NULL);
AUTO_ASSERT(0 == shape->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}

void np02_loc_grid::get_shapes_near_shape(const np02_shape *s,
    shape_vec *shapes) const{
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
np02_shape *shape;
AUTO_ASSERT(m_idx_shape_vec.empty());
shp_type_idx_shape_vec *idx_shape_vec =
    const_cast<shp_type_idx_shape_vec *>(&m_idx_shape_vec);
const np02_loc_grid_node *lg_node =
    (NULL == s) ? NULL : s->get_head_loc_grid_node();
if((NULL != s) && (NULL != shapes) && (NULL != lg_node)){
    AUTO_ASSERT(this == s->get_loc_grid());

    /* search all grid squares occupied by s */
    while(NULL != lg_node){
        AUTO_ASSERT(0 == lg_node->verify_data(AA_ERR_BUF(),
            AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR()));

        /* search up in same grid square */
        const np02_loc_grid_node *lg_sq_node = lg_node->get_prev();
        while(NULL != lg_sq_node){
            shape = lg_sq_node->get_owner();
            AUTO_ASSERT(shape != s);
            idx_shape_vec->push_back(shp_type_idx_shape(shp_type_idx_pair(
                shape->get_shape_type(),shape->get_shape_idx()), shape) );
            AA_XDBG_ASSERT(0 == lg_sq_node->verify_data(AA_ERR_BUF(),
                AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR()),
                CF01_AA_DEBUG_LEVEL_2);
            lg_sq_node = lg_sq_node->get_prev();
            }

        /* search down in same grid square */
        lg_sq_node = lg_node->get_next();
        while(NULL != lg_sq_node){
            shape = lg_sq_node->get_owner();
            AUTO_ASSERT(shape != s);
            idx_shape_vec->push_back(shp_type_idx_shape(shp_type_idx_pair(
                shape->get_shape_type(),shape->get_shape_idx()), shape) );
            AA_XDBG_ASSERT(0 == lg_sq_node->verify_data(AA_ERR_BUF(),
                AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR()),
                CF01_AA_DEBUG_LEVEL_2);
            lg_sq_node = lg_sq_node->get_next();
            }

        lg_node = lg_node->get_s_next();
        }

    /* sort search vec, remove duplicates */
    std::sort(idx_shape_vec->begin(), idx_shape_vec->end());
    shp_type_idx_shape_vec_itr unique_end_itr =
        std::unique(idx_shape_vec->begin(), idx_shape_vec->end());

    /* copy nearby shapes to output vector */
    shp_type_idx_shape_vec_itr unique_itr = idx_shape_vec->begin();
    for(; unique_itr != unique_end_itr; ++unique_itr){
        shape = unique_itr->second;
        shapes->push_back(shape);
        }

    /* restore utility vector */
    idx_shape_vec->clear();
    }
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}

void np02_loc_grid::get_shapes_in_bb(const np02_xy& xy_min,
    const np02_xy& xy_max, shape_vec *shapes ) const{
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_ALWAYS_ASSERT(xy_min.get_x() <= xy_max.get_x());
AA_ALWAYS_ASSERT(xy_min.get_y() <= xy_max.get_y());
AA_ALWAYS_ASSERT(NULL != shapes);
np02_shape *shape;

/* utility vector */
AUTO_ASSERT(m_idx_shape_vec.empty());
shp_type_idx_shape_vec *idx_shape_vec =
    const_cast<shp_type_idx_shape_vec *>(&m_idx_shape_vec);

/* search locator grid, fill search_vec */
size_t lg_idx;
uint16_t i, j;
np02_uint16_pair ij_min, ij_max;
m_loc_grid_dim.get_bb_indices(xy_min, xy_max, &ij_min, &ij_max);
for(i = ij_min.first; i <= ij_max.first; ++i){
    for(j = ij_min.second; j <= ij_max.second; ++j){
        lg_idx = static_cast<size_t>(j) + (static_cast<size_t>(i) *
            static_cast<size_t>(m_loc_grid_dim.get_h()));
        AA_ALWAYS_ASSERT(lg_idx < m_loc_grid_vec.size());
        if(lg_idx < m_loc_grid_vec.size()){
            np02_loc_grid_node *lg_node = m_loc_grid_vec[lg_idx];
            while(NULL != lg_node){
                shape = lg_node->get_owner();
                AA_ALWAYS_ASSERT(NULL != shape);
                idx_shape_vec->push_back(shp_type_idx_shape(shp_type_idx_pair(
                    shape->get_shape_type(),shape->get_shape_idx()), shape) );
                lg_node = lg_node->get_next();
                }
            }
        }
    }

/* sort search vec, remove duplicates */
std::sort(idx_shape_vec->begin(), idx_shape_vec->end());
shp_type_idx_shape_vec_itr unique_end_itr =
    std::unique(idx_shape_vec->begin(), idx_shape_vec->end());

/* copy shapes to output vector */
shp_type_idx_shape_vec_itr unique_itr = idx_shape_vec->begin();
for(; unique_itr != unique_end_itr; ++unique_itr){
    shape = unique_itr->second;
    shapes->push_back(shape);
    }

/* restore utility vector */
idx_shape_vec->clear();

AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}

int np02_loc_grid::verify_data( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;

size_t count, loc_grid_node_vec_sz_check;
uint16_t w=0, h=0, i=0, j=0;
size_t gn_vec_idx=0;
np02_shape *shape=NULL;
np02_loc_grid_node *loc_grid_node=NULL;

w=m_loc_grid_dim.get_w();
h=m_loc_grid_dim.get_h();
err_cnt=m_loc_grid_dim.verify_data(err_msg,err_msg_capacity,err_msg_pos);
loc_grid_node_vec_sz_check = static_cast<size_t>(w) * static_cast<size_t>(h);
size_t max_expected_shape_count= 1000000;
const np02_shp_alloc *shp_alloc = get_shp_alloc();
if(NULL != shp_alloc){
    max_expected_shape_count = shp_alloc->alloc_get_circle_count() +
        shp_alloc->alloc_get_line_seg_count() +
        shp_alloc->alloc_get_rect_count() +
        shp_alloc->alloc_get_polygon_count();

    if( this != shp_alloc->alloc_get_loc_grid_by_idx(m_alloc_idx)){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "loc_grid=%x  !=  shp_alloc=%x->get loc_grid by "
            "(idx=%i) = %x\n",
            this, shp_alloc, m_alloc_idx,
            shp_alloc->alloc_get_loc_grid_by_idx(m_alloc_idx) );
        }
    }

if(m_loc_grid_vec.size() != loc_grid_node_vec_sz_check){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "loc_grid=%x m_loc_grid_vec.size=%i != check_size=%i\n",
        this, m_loc_grid_vec.size(), loc_grid_node_vec_sz_check);
    }

if( m_extra_search_d < 0.0 ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "loc_grid=%x extra_search_d=%f\n", this, m_extra_search_d);
    }

for(gn_vec_idx = 0; gn_vec_idx<m_loc_grid_vec.size(); ++gn_vec_idx){
    i = (h > 0) ? static_cast<uint16_t>(gn_vec_idx / h) : 0;
    if(i > w){i = (w>0) ? (w-1) : 0; }
    j = (h > 0) ? static_cast<uint16_t>(gn_vec_idx % h) : 0;

    /* add every loc grid node to total_gn_vec_g */
    loc_grid_node = m_loc_grid_vec.at(gn_vec_idx);
    count = 0;
    while((NULL != loc_grid_node) && (count <= max_expected_shape_count)){
        ++count;
        err_cnt += loc_grid_node->verify_data(err_msg, err_msg_capacity,
            err_msg_pos);
        if(this != loc_grid_node->get_loc_grid()){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "loc_grid=%x  this != (loc_grid_node=%x:i=%i:j=%i)"
                "->loc_grid=%x\n",
                this, loc_grid_node, i, j, loc_grid_node->get_loc_grid() ); }
        if(loc_grid_node->get_i() != i){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "loc_grid=%x (loc_grid_node=%x)->i=%i != i=%i\n",
                this, loc_grid_node, loc_grid_node->get_i(), i ); }
        if(loc_grid_node->get_j() != j){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "loc_grid=%x (loc_grid_node=%x)->j=%i != j=%i\n",
                this, loc_grid_node, loc_grid_node->get_j(), j ); }
        loc_grid_node = loc_grid_node->get_next();
        }
    if(count > max_expected_shape_count){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "loc_grid=%x i=%i j=%i, count=%i > max_expected_count=%i\n",
            this, i, j, count, max_expected_shape_count ); 
        }
    }

return err_cnt;
}

std::ostream& np02_loc_grid::ostream_output(std::ostream& os) const{
os << "<loc_grid=" << std::hex << this << std::dec << ">\n";
os << "<alloc_idx>" << m_alloc_idx << "</alloc_idx>\n";
//os << "<owner>" << std::hex << m_owner << std::dec << "</owner>\n";
os << "<extra_search_d>" << m_extra_search_d << "</extra_search_d>\n";
m_loc_grid_dim.ostream_output(os); 
os << "<count>" << m_loc_grid_vec.size() << "</count>\n";
uint16_t i = 0, j = 0;
for(loc_grid_node_vec_citr lg_node_itr=m_loc_grid_vec.begin();
    lg_node_itr!=m_loc_grid_vec.end(); ++lg_node_itr){
    const np02_loc_grid_node *lg_node = *lg_node_itr;
    os << "</i=" << i << "></j=" << j << ">";
    if(NULL == lg_node){
        os << "</lg_node=NULL>\n";
        }
    else{
        const size_t lg_count_max = 1000000;
        size_t lg_count = 0;
        while((lg_count < lg_count_max) && (NULL != lg_node)){
            ++lg_count;
            lg_node->ostream_output(os);
            lg_node = lg_node->get_next();
            }
        }
    ++j;
    if (j >= m_loc_grid_dim.get_h()){ j = 0; ++i; }
    }
os << "</loc_grid>\n";
return os;
}

std::ostream& np02_loc_grid::ostream_output_brief_table(
    std::ostream& os) const{
os << "<loc_grid=" << std::hex << this << std::dec << ">\n";
os << "<alloc_idx>" << m_alloc_idx << "</alloc_idx>\n";
//os << "<owner>" << std::hex << m_owner << std::dec << "</owner>\n";
os << "<extra_search_d>" << m_extra_search_d << "</extra_search_d>\n";
m_loc_grid_dim.ostream_output(os); 
os << "<count>" << m_loc_grid_vec.size() << "</count>\n";
uint16_t i = 0, j = 0;
for(loc_grid_node_vec_citr lg_node_itr=m_loc_grid_vec.begin();
    lg_node_itr!=m_loc_grid_vec.end(); ++lg_node_itr){
    const np02_loc_grid_node *lg_node = *lg_node_itr;
    os << "</i=" << i << "></j=" << j << ">";
    if(NULL == lg_node){
        os << "</lg_node=NULL>\n";
        }
    else{
        const size_t lg_count_max = 1000000;
        size_t lg_count = 0;
        while((lg_count < lg_count_max) && (NULL != lg_node)){
            ++lg_count;
            lg_node = lg_node->get_next();
            }
        os << "<node_count>" << lg_count << "<node_count>\n";
        }
    ++j;
    if (j >= m_loc_grid_dim.get_h()){ j = 0; ++i; }
    }
os << "</loc_grid>\n";
return os;
}


void np02_loc_grid::write_bmp_file(const np02_xy& xy_min,
    const double& pixel_num, const np02_bmp_color& color,
    np02_bmp_file *bmp_file) const{
write_bmp_file_grid_lines(xy_min, pixel_num, color, bmp_file);
write_bmp_file_grid_info(xy_min, pixel_num, color, bmp_file);
}

void np02_loc_grid::write_bmp_file_grid_lines(const np02_xy& xy_min,
    const double& pixel_num, const np02_bmp_color& color,
    np02_bmp_file *bmp_file) const{
m_loc_grid_dim.write_bmp_file(xy_min, pixel_num, color, bmp_file);
}

void np02_loc_grid::write_bmp_file_grid_info(const np02_xy& xy_min,
    const double& pixel_num, const np02_bmp_color& color,
    np02_bmp_file *bmp_file) const{
if(NULL != bmp_file){
    const double hh = (m_loc_grid_dim.get_sq_size() * pixel_num) / 12;
    double ii0,ii1,jj0,jj1;
    uint16_t i = 0, j = 0;
    for(loc_grid_node_vec_citr lg_node_itr=m_loc_grid_vec.begin();
        lg_node_itr!=m_loc_grid_vec.end(); ++lg_node_itr){
        const np02_loc_grid_node *lg_node = *lg_node_itr;
        size_t shape_count = 0;
        if(NULL != lg_node){
            const size_t shape_count_max = 1000000;
            while((shape_count < shape_count_max) && (NULL != lg_node)){
                ++shape_count;
                lg_node = lg_node->get_next();
                }
            }

        char msg_buf[32];
        sprintf(msg_buf, "i%i j%i s%i", i,j,shape_count);

        ii0 = (pixel_num * (m_loc_grid_dim.get_x_min() + 
            (static_cast<double>(i)* m_loc_grid_dim.get_sq_size()) - 
            xy_min.get_x())) + (hh/2.0);
        jj0 = (pixel_num * (m_loc_grid_dim.get_y_min() +
            (static_cast<double>(j)* m_loc_grid_dim.get_sq_size()) -
            xy_min.get_x())) + (hh/2.0);

        bmp_file->draw_text(msg_buf, ii0, jj0, hh, 0.0, color);

        ++j;
        if (j >= m_loc_grid_dim.get_h()){ j = 0; ++i; }
        }

    /* extra search distance bar */
    ii0 = (pixel_num * (m_loc_grid_dim.get_x_min() -  xy_min.get_x()));
    jj0 = (pixel_num * (m_loc_grid_dim.get_y_min() - xy_min.get_x())) - 
        (7.0*hh/8.0);
    ii1 = (pixel_num * (m_loc_grid_dim.get_x_min() + m_extra_search_d - 
        xy_min.get_x()));
    jj1 = (pixel_num * (m_loc_grid_dim.get_y_min() - xy_min.get_x()))-(hh/8.0);
    bmp_file->draw_box( static_cast<int32_t>(ii0), static_cast<int32_t>(jj0),
        static_cast<int32_t>(ii1), static_cast<int32_t>(jj1), color);

    ii0 = (pixel_num * (m_loc_grid_dim.get_x_min() + m_extra_search_d - 
        xy_min.get_x()));
    jj0 = (pixel_num * (m_loc_grid_dim.get_y_min() - xy_min.get_x())) - hh;
    bmp_file->draw_text("extra_search_d", ii0, jj0, hh, 0.0, color);
    }
}


np02_loc_grid_node *np02_loc_grid::alloc_loc_grid_node(){
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
np02_loc_grid_node *n = NULL;
np02_shp_alloc *w = get_shp_alloc();
if( NULL == w ){ n = new np02_loc_grid_node();  }
else{ n = w->alloc_loc_grid_node(); }
if(NULL != n){
    n->set_loc_grid(this);
    AUTO_ASSERT(NULL == n->get_owner());
    AUTO_ASSERT(0 == n->get_i());
    AUTO_ASSERT(0 == n->get_j());
    AUTO_ASSERT(NULL == n->get_prev());
    AUTO_ASSERT(NULL == n->get_next());
    AUTO_ASSERT(NULL == n->get_s_prev());
    AUTO_ASSERT(NULL == n->get_s_next());
    AUTO_ASSERT(0 == n->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
    }
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
return n;
}

void np02_loc_grid::free_loc_grid_node(np02_loc_grid_node *n){
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != n);
AUTO_ASSERT(this == n->get_loc_grid());

np02_shp_alloc *w = get_shp_alloc();
if(NULL == w){
    delete n;
    }
else{
    n->set_owner(NULL);
    n->set_loc_grid(NULL);
    n->set_i(0);
    n->set_j(0);
    n->set_prev(NULL);
    n->set_next(NULL);
    n->set_s_prev(NULL);
    n->set_s_next(NULL);
    w->free_loc_grid_node(n);
    AUTO_ASSERT(0 == n->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
        AA_ERR_BUF_POS_PTR() ));
    }
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}

void np02_loc_grid::clear_loc_grid(){
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
loc_grid_node_vec_citr lg_itr = m_loc_grid_vec.begin();
for(; lg_itr != m_loc_grid_vec.end(); ++lg_itr){
    np02_loc_grid_node *lg_node = *lg_itr;
    if(NULL != lg_node){
        static const size_t max_shape_count = 10000000;
        size_t shape_count = 0;
        while((NULL != lg_node) && (shape_count < max_shape_count)){
            ++shape_count;
            np02_shape *shape = lg_node->get_owner();
            AA_ALWAYS_ASSERT(NULL != shape);
            remove_shape_from_loc_grid(shape);
            lg_node = *lg_itr;
            }
        AUTO_ASSERT(shape_count < max_shape_count);
        }
    }
m_loc_grid_vec.clear();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}











np02_shp_alloc::np02_shp_alloc():
    /* allocator */
    m_alloc_circle_vec(),
    m_circle_free_chain(NULL),
    m_alloc_line_seg_vec(),
    m_line_seg_free_chain(NULL),
    m_alloc_rect_vec(),
    m_rect_free_chain(NULL),
    m_alloc_polygon_vec(),
    m_polygon_free_chain(NULL),
    m_alloc_loc_grid_node_vec(),
    m_loc_grid_node_free_chain(NULL),
    m_alloc_loc_grid_vec()
{
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()));
AA_DECR_CALL_DEPTH();
}

np02_shp_alloc::~np02_shp_alloc(){
alloc_delete_all();
}


void np02_shp_alloc::free_shape(np02_shape *shape)
{
AA_INCR_CALL_DEPTH();
if(NULL != shape){
    np02_circle *circle;
    np02_line_seg *line_seg;
    np02_rect *rect;
    np02_polygon *polygon;
    if( (circle = dynamic_cast<np02_circle *>(shape)) != NULL ){
        free_circle(circle);
        }
    else if((line_seg=dynamic_cast<np02_line_seg *>(shape)) != NULL){
        free_line_seg(line_seg);
        }
    else if((rect=dynamic_cast<np02_rect *>(shape)) != NULL){
        free_rect(rect);
        }
    else if((polygon=dynamic_cast<np02_polygon *>(polygon)) != NULL){
        free_polygon(polygon);
        }
    else{ AA_ALWAYS_ASSERT(false); }
    }
AA_DECR_CALL_DEPTH();
}

np02_circle *np02_shp_alloc::alloc_circle(){
AA_INCR_CALL_DEPTH();
np02_circle *circle = NULL;
if(NULL == m_circle_free_chain){
    circle = new np02_circle();
    circle->set_shape_idx(m_alloc_circle_vec.size());
    circle->set_shp_alloc(this);
    m_alloc_circle_vec.push_back(circle);
    }
else{
    circle = m_circle_free_chain;
    m_circle_free_chain = m_circle_free_chain->get_free_chain_next();
    circle->set_free_chain_next(NULL);
    }
AUTO_ASSERT(0 == circle->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_DECR_CALL_DEPTH();
return circle;
}

void np02_shp_alloc::free_circle(np02_circle *circle){
AA_ALWAYS_ASSERT( NULL != circle );
AA_ALWAYS_ASSERT( this == circle->get_shp_alloc() );
circle->destruct();
if(NULL != m_circle_free_chain){
    circle->set_free_chain_next(m_circle_free_chain);
    }
m_circle_free_chain = circle;
circle->set_shp_owner_idx(NP02_SHP_OWNER_INVALID_IDX);
}

np02_line_seg *np02_shp_alloc::alloc_line_seg(){
AA_INCR_CALL_DEPTH();
np02_line_seg *line_seg = NULL;
if(NULL == m_line_seg_free_chain){
    line_seg = new np02_line_seg();
    line_seg->set_shape_idx(m_alloc_line_seg_vec.size());
    line_seg->set_shp_alloc(this);
    m_alloc_line_seg_vec.push_back(line_seg);
    }
else{
    line_seg = m_line_seg_free_chain;
    m_line_seg_free_chain = m_line_seg_free_chain->get_free_chain_next();
    line_seg->set_free_chain_next(NULL);
    }
AUTO_ASSERT(0 == line_seg->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_DECR_CALL_DEPTH();
return line_seg;
}

void np02_shp_alloc::free_line_seg(np02_line_seg *line_seg){
AA_ALWAYS_ASSERT( NULL != line_seg );
AA_ALWAYS_ASSERT( this == line_seg->get_shp_alloc() );
line_seg->destruct();
if(NULL != m_line_seg_free_chain){
    line_seg->set_free_chain_next(m_line_seg_free_chain);
    }
m_line_seg_free_chain = line_seg;
line_seg->set_shp_owner_idx(NP02_SHP_OWNER_INVALID_IDX);
}

np02_rect *np02_shp_alloc::alloc_rect(){
AA_INCR_CALL_DEPTH();
np02_rect *rect = NULL;
if(NULL == m_rect_free_chain){
    rect = new np02_rect();
    rect->set_shape_idx(m_alloc_rect_vec.size());
    rect->set_shp_alloc(this);
    m_alloc_rect_vec.push_back(rect);
    }
else{
    rect = m_rect_free_chain;
    m_rect_free_chain = m_rect_free_chain->get_free_chain_next();
    rect->set_free_chain_next(NULL);
    }
AUTO_ASSERT(0 == rect->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_DECR_CALL_DEPTH();
return rect;
}

void np02_shp_alloc::free_rect(np02_rect *rect){
AA_ALWAYS_ASSERT( NULL != rect );
AA_ALWAYS_ASSERT( this == rect->get_shp_alloc() );
rect->destruct();
if(NULL != m_rect_free_chain){
    rect->set_free_chain_next(m_rect_free_chain);
    }
m_rect_free_chain = rect;
rect->set_shp_owner_idx(NP02_SHP_OWNER_INVALID_IDX);
}

np02_polygon *np02_shp_alloc::alloc_polygon(){
AA_INCR_CALL_DEPTH();
np02_polygon *polygon = NULL;
if(NULL == m_polygon_free_chain){
    polygon = new np02_polygon();
    polygon->set_shape_idx(m_alloc_polygon_vec.size());
    polygon->set_shp_alloc(this);
    m_alloc_polygon_vec.push_back(polygon);
    }
else{
    polygon = m_polygon_free_chain;
    m_polygon_free_chain = m_polygon_free_chain->get_free_chain_next();
    polygon->set_free_chain_next(NULL);
    }
AUTO_ASSERT(0 == polygon->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_DECR_CALL_DEPTH();
return polygon;
}

void np02_shp_alloc::free_polygon(np02_polygon *polygon){
AA_ALWAYS_ASSERT( NULL != polygon );
AA_ALWAYS_ASSERT( this == polygon->get_shp_alloc() );
polygon->destruct();
if(NULL != m_polygon_free_chain){
    polygon->set_free_chain_next(m_polygon_free_chain);
    }
m_polygon_free_chain = polygon;
polygon->set_shp_owner_idx(NP02_SHP_OWNER_INVALID_IDX);
}

np02_loc_grid_node *np02_shp_alloc::alloc_loc_grid_node(){
AA_INCR_CALL_DEPTH();
np02_loc_grid_node *node = NULL;
if(NULL == m_loc_grid_node_free_chain){
    node = new np02_loc_grid_node();
    node->set_alloc_idx(m_alloc_loc_grid_node_vec.size());
    m_alloc_loc_grid_node_vec.push_back(node);
    }
else{
    node = m_loc_grid_node_free_chain;
    m_loc_grid_node_free_chain = m_loc_grid_node_free_chain->get_free_chain_next();
    node->set_free_chain_next(NULL);
    }
AUTO_ASSERT(0 == node->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_DECR_CALL_DEPTH();
return node;
}

void np02_shp_alloc::free_loc_grid_node(np02_loc_grid_node *loc_grid_node){
/* TODO: call node->destruct() (might be named differently.) to free resources */
if(NULL != m_loc_grid_node_free_chain){
    loc_grid_node->set_free_chain_next(m_loc_grid_node_free_chain);
    }
m_loc_grid_node_free_chain = loc_grid_node;
loc_grid_node->set_owner(NULL);
}

np02_loc_grid *np02_shp_alloc::alloc_loc_grid(){
AA_INCR_CALL_DEPTH();
np02_loc_grid *loc_grid = new np02_loc_grid();
loc_grid->set_shp_alloc(this);
loc_grid->set_alloc_idx(m_alloc_loc_grid_vec.size());
m_alloc_loc_grid_vec.push_back(loc_grid);
AUTO_ASSERT(0 == loc_grid->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_DECR_CALL_DEPTH();
return loc_grid;
}

void np02_shp_alloc::free_loc_grid(np02_loc_grid *loc_grid){
if(NULL != loc_grid){
    const size_t& alloc_idx = loc_grid->get_alloc_idx();
    if(alloc_idx >= m_alloc_loc_grid_vec.size()){
        AA_ALWAYS_ASSERT(false);
        }
    else{
        AUTO_ASSERT(m_alloc_loc_grid_vec.at(alloc_idx) == loc_grid);
        m_alloc_loc_grid_vec.at(alloc_idx)=NULL;
        }
    delete loc_grid;
    }
}

int np02_shp_alloc::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const{
size_t i;
int err_cnt = 0;
err_cnt += verify_data_alloc_circle(err_msg, err_msg_capacity, err_msg_pos );
err_cnt += verify_data_alloc_line_seg(err_msg, err_msg_capacity, err_msg_pos);
err_cnt += verify_data_alloc_rect(err_msg, err_msg_capacity, err_msg_pos );
err_cnt += verify_data_alloc_polygon(err_msg, err_msg_capacity, err_msg_pos );
err_cnt+=verify_data_alloc_loc_grid_node(err_msg,err_msg_capacity,err_msg_pos);
err_cnt += verify_data_alloc_loc_grid(err_msg, err_msg_capacity, err_msg_pos );
return err_cnt;
}

std::ostream& np02_shp_alloc::ostream_output(std::ostream& os) const{
os << "<shp_alloc><this>" << std::hex << this << std::dec << "</this>\n";

/* allocator */
os << "</alloc_circle_count=" << m_alloc_circle_vec.size() << ">" 
    //<< "</circle_free_chain_size=" << get_circle_free_chain_size() << ">"
    << "\n";
os << "</alloc_line_seg_count=" << m_alloc_line_seg_vec.size() << ">" 
    //<< "</cline_seg_free_chain_size=" << get_line_seg_free_chain_size() << ">"
    << "\n";
os << "</alloc_rect_count=" << m_alloc_rect_vec.size() << ">" 
    //<< "</rect_free_chain_size=" << get_rect_free_chain_size() << ">"
    << "\n";
os << "</alloc_polygon_count=" << m_alloc_polygon_vec.size() << ">" 
    //<< "</polygon_free_chain_size=" << get_polygon_free_chain_size() << ">"
    << "\n";
os << "</alloc_loc_grid_node_count=" << m_alloc_loc_grid_node_vec.size() << ">" 
    //<< "</loc_grid_node_free_chain_size=" << get_loc_grid_node_free_chain_size() << ">"
    << "\n";
os << "</alloc_loc_grid_count=" << m_alloc_loc_grid_vec.size() << ">" << "\n";
os << "</shp_alloc>\n";

return os;
}

/* destructor implementation.  Delete all objects. */
void np02_shp_alloc::alloc_delete_all(){
AA_INCR_CALL_DEPTH();
while (NULL != m_circle_free_chain) {
    np02_circle *circle = m_circle_free_chain;
    m_circle_free_chain = m_circle_free_chain->get_free_chain_next();
    circle->set_free_chain_next(NULL);
    delete circle;
    }
m_alloc_circle_vec.clear();

m_line_seg_free_chain = NULL;
while (NULL != m_line_seg_free_chain) {
    np02_line_seg *line_seg = m_line_seg_free_chain;
    m_line_seg_free_chain = m_line_seg_free_chain->get_free_chain_next();
    line_seg->set_free_chain_next(NULL);
    delete line_seg;
    }
m_alloc_line_seg_vec.clear();

while (NULL != m_rect_free_chain) {
    np02_rect *rect = m_rect_free_chain;
    m_rect_free_chain = m_rect_free_chain->get_free_chain_next();
    rect->set_free_chain_next(NULL);
    delete rect;
    }
m_alloc_rect_vec.clear();

while (NULL != m_polygon_free_chain) {
    np02_polygon *polygon = m_polygon_free_chain;
    m_polygon_free_chain = m_polygon_free_chain->get_free_chain_next();
    polygon->set_free_chain_next(NULL);
    delete polygon;
    }
m_alloc_polygon_vec.clear();

while (NULL != m_loc_grid_node_free_chain) {
    np02_loc_grid_node *loc_grid_node = m_loc_grid_node_free_chain;
    m_loc_grid_node_free_chain =
        m_loc_grid_node_free_chain->get_free_chain_next();
    loc_grid_node->set_free_chain_next(NULL);
    delete loc_grid_node;
    }
m_alloc_loc_grid_node_vec.clear();

for(loc_grid_vec::const_iterator loc_grid_itr = m_alloc_loc_grid_vec.begin();
    loc_grid_itr != m_alloc_loc_grid_vec.end(); ++loc_grid_itr){
    if(NULL != *loc_grid_itr){
        delete *loc_grid_itr;
        }
    }
m_alloc_loc_grid_vec.clear();

AA_DECR_CALL_DEPTH();
}

int np02_shp_alloc::verify_data_alloc_circle( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;
size_t i;
const np02_circle *circle = NULL;
for(i = 0; i < m_alloc_circle_vec.size(); ++i){
    circle = m_alloc_circle_vec.at(i);
    if(NULL == circle){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  m_alloc_circle_vec.at(%i)=NULL\n",
            this, i ); 
        }
    else{
        if(circle->get_shape_idx() != i){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  circle(%i)=%x->shape_idx=%i\n",
                this, i, circle, circle->get_shape_idx() );
            }
        const np02_shp_alloc *circle_shp_alloc = circle->get_shp_alloc();
        if((NULL != circle_shp_alloc) && (this != circle_shp_alloc)){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  != circle(%i)=%x -> shp_alloc=%x\n",
                this, i, circle, circle_shp_alloc ); 
            }
        }
    }

i=0;
circle=m_circle_free_chain;
while((NULL != circle) && (i < m_alloc_circle_vec.size())){
    const size_t circle_shp_idx = circle->get_shape_idx();
    if(m_alloc_circle_vec.size() <= circle_shp_idx){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain circle %i "
            "shp_idx=%i >= vec sz=%i\n", this, i, circle_shp_idx,
            m_alloc_circle_vec.size() ); 
        }
    else if(m_alloc_circle_vec.at(circle_shp_idx) != circle){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain circle %i =%x"
            " != vec.at(shp_idx=%i) = %x\n", this, i, circle, circle_shp_idx,
            m_alloc_circle_vec.at(circle_shp_idx) ); 
        }
    
    ++i;
    circle = circle->get_free_chain_next();
    }
return err_cnt;
}

int np02_shp_alloc::verify_data_alloc_line_seg( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;
size_t i;
const np02_line_seg *line_seg = NULL;
for(i = 0; i < m_alloc_line_seg_vec.size(); ++i){
    line_seg = m_alloc_line_seg_vec.at(i);
    if(NULL == line_seg){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  m_alloc_line_seg_vec.at(%i)=NULL\n",
            this, i ); 
        }
    else{
        if(line_seg->get_shape_idx() != i){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  line_seg(%i)=%x->shape_idx=%i\n",
                this, i, line_seg, line_seg->get_shape_idx() );
            }
        const np02_shp_alloc *line_seg_shp_alloc = line_seg->get_shp_alloc();
        if((NULL != line_seg_shp_alloc) && (this != line_seg_shp_alloc)){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  != line_seg(%i)=%x -> shp_alloc=%x\n",
                this, i, line_seg, line_seg_shp_alloc ); 
            }
        }
    }

i=0;
line_seg=m_line_seg_free_chain;
while((NULL != line_seg) && (i < m_alloc_line_seg_vec.size())){
    const size_t line_seg_shp_idx = line_seg->get_shape_idx();
    if(m_alloc_line_seg_vec.size() <= line_seg_shp_idx){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain line_seg %i "
            "shp_idx=%i >= vec sz=%i\n", this, i, line_seg_shp_idx,
            m_alloc_line_seg_vec.size() ); 
        }
    else if(m_alloc_line_seg_vec.at(line_seg_shp_idx) != line_seg){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain line_seg %i =%x"
            " != vec.at(shp_idx=%i) = %x\n", this, i, line_seg, line_seg_shp_idx,
            m_alloc_line_seg_vec.at(line_seg_shp_idx) ); 
        }
    
    ++i;
    line_seg = line_seg->get_free_chain_next();
    }
return err_cnt;
}

int np02_shp_alloc::verify_data_alloc_rect( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;
size_t i;
const np02_rect *rect = NULL;
for(i = 0; i < m_alloc_rect_vec.size(); ++i){
    rect = m_alloc_rect_vec.at(i);
    if(NULL == rect){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  m_alloc_rect_vec.at(%i)=NULL\n",
            this, i ); 
        }
    else{
        if(rect->get_shape_idx() != i){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  rect(%i)=%x->shape_idx=%i\n",
                this, i, rect, rect->get_shape_idx() );
            }
        const np02_shp_alloc *rect_shp_alloc = rect->get_shp_alloc();
        if((NULL != rect_shp_alloc) && (this != rect_shp_alloc)){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  != rect(%i)=%x -> shp_alloc=%x\n",
                this, i, rect, rect_shp_alloc ); 
            }
        }
    }

i=0;
rect=m_rect_free_chain;
while((NULL != rect) && (i < m_alloc_rect_vec.size())){
    const size_t rect_shp_idx = rect->get_shape_idx();
    if(m_alloc_rect_vec.size() <= rect_shp_idx){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain rect %i "
            "shp_idx=%i >= vec sz=%i\n", this, i, rect_shp_idx,
            m_alloc_rect_vec.size() ); 
        }
    else if(m_alloc_rect_vec.at(rect_shp_idx) != rect){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain rect %i =%x"
            " != vec.at(shp_idx=%i) = %x\n", this, i, rect, rect_shp_idx,
            m_alloc_rect_vec.at(rect_shp_idx) ); 
        }
    
    ++i;
    rect = rect->get_free_chain_next();
    }
return err_cnt;
}

int np02_shp_alloc::verify_data_alloc_polygon( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;
size_t i;
const np02_polygon *polygon = NULL;
for(i = 0; i < m_alloc_polygon_vec.size(); ++i){
    polygon = m_alloc_polygon_vec.at(i);
    if(NULL == polygon){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  m_alloc_polygon_vec.at(%i)=NULL\n",
            this, i ); 
        }
    else{
        if(polygon->get_shape_idx() != i){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  polygon(%i)=%x->shape_idx=%i\n",
                this, i, polygon, polygon->get_shape_idx() );
            }
        const np02_shp_alloc *polygon_shp_alloc = polygon->get_shp_alloc();
        if((NULL != polygon_shp_alloc) && (this != polygon_shp_alloc)){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  != polygon(%i)=%x -> shp_alloc=%x\n",
                this, i, polygon, polygon_shp_alloc ); 
            }
        }
    }

i=0;
polygon=m_polygon_free_chain;
while((NULL != polygon) && (i < m_alloc_polygon_vec.size())){
    const size_t polygon_shp_idx = polygon->get_shape_idx();
    if(m_alloc_polygon_vec.size() <= polygon_shp_idx){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain polygon %i "
            "shp_idx=%i >= vec sz=%i\n", this, i, polygon_shp_idx,
            m_alloc_polygon_vec.size() ); 
        }
    else if(m_alloc_polygon_vec.at(polygon_shp_idx) != polygon){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain polygon %i =%x"
            " != vec.at(shp_idx=%i) = %x\n", this, i, polygon, polygon_shp_idx,
            m_alloc_polygon_vec.at(polygon_shp_idx) ); 
        }
    
    ++i;
    polygon = polygon->get_free_chain_next();
    }
return err_cnt;
}

int np02_shp_alloc::verify_data_alloc_loc_grid_node( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;
size_t i;
const np02_loc_grid_node *lg_node = NULL;
for(i = 0; i < m_alloc_loc_grid_node_vec.size(); ++i){
    lg_node = m_alloc_loc_grid_node_vec.at(i);
    if(NULL == lg_node){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  m_alloc_loc_grid_node_vec.at(%i)=NULL\n",
            this, i ); 
        }
    else{
        if(lg_node->get_alloc_idx() != i){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  lg_node(%i)=%x->alloc_idx=%i\n",
                this, i, lg_node, lg_node->get_alloc_idx() );
            }
        const np02_shp_alloc *lg_node_shp_alloc = lg_node->get_shp_alloc();
        if((NULL != lg_node_shp_alloc) && (this != lg_node_shp_alloc)){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  != lg_node(%i)=%x -> shp_alloc=%x\n",
                this, i, lg_node, lg_node_shp_alloc ); 
            }
        }
    }

i=0;
lg_node=m_loc_grid_node_free_chain;
while((NULL != lg_node) && (i < m_alloc_loc_grid_node_vec.size())){
    const size_t alloc_idx = lg_node->get_alloc_idx();
    if(m_alloc_loc_grid_node_vec.size() <= alloc_idx){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain lg_node %i "
            "alloc_idx=%i >= vec sz=%i\n", this, i, alloc_idx,
            m_alloc_loc_grid_node_vec.size() ); 
        }
    else if(m_alloc_loc_grid_node_vec.at(alloc_idx) != lg_node){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain lg_node %i =%x"
            " != vec.at(alloc_idx=%i) = %x\n", this, i, lg_node, alloc_idx,
            m_alloc_loc_grid_node_vec.at(alloc_idx) ); 
        }
    
    ++i;
    lg_node = lg_node->get_free_chain_next();
    }
return err_cnt;
}

int np02_shp_alloc::verify_data_alloc_loc_grid( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;
size_t i;
const np02_loc_grid *loc_grid = NULL;
for(i = 0; i < m_alloc_loc_grid_vec.size(); ++i){
    loc_grid = m_alloc_loc_grid_vec.at(i);
    if(NULL != loc_grid){
        if(loc_grid->get_alloc_idx() != i){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  loc_grid(%i)=%x->alloc_idx=%i\n",
                this, i, loc_grid, loc_grid->get_alloc_idx() );
            }
        const np02_shp_alloc *loc_grid_shp_alloc = loc_grid->get_shp_alloc();
        if(this != loc_grid_shp_alloc){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  != loc_grid(%i)=%x -> shp_alloc=%x\n",
                this, i, loc_grid, loc_grid_shp_alloc ); 
            }
        }
    }

return err_cnt;
}

int np02_shape_test::run_shape_test(const int& shape_test_number,const int&
    shape_test_iteration_count, const int& shape_test_rand_seed){
np02_shape_test t;
t.set_shape_test_number( shape_test_number );
t.set_shape_test_iteration_count(shape_test_iteration_count);
t.set_shape_test_rand_seed( shape_test_rand_seed );
const int error_code = t.execute();
return error_code;
}

np02_shape_test::np02_shape_test(): m_shape_test_number(1),
    m_shape_test_iteration_count(1), m_shape_test_rand_seed(1),
    m_iteration(0), m_rand_uint32(4294967291ul)
{}

np02_shape_test::~np02_shape_test(){}

void np02_shape_test::advance_rand(){
/* advance random number, pg 359, The Standard C Library (c)1991,
P. J. Plauger, ISBN 0-13-131509-9 */
m_rand_uint32= (m_rand_uint32 * 1103515245) + 12345;

/* 
https://en.wikipedia.org/wiki/Linear-feedback_shift_register#Xorshift_LFSRs 
*/
m_rand_uint32 ^= ( m_rand_uint32 >> 7 );
m_rand_uint32 ^= ( m_rand_uint32 << 9 );
m_rand_uint32 ^= ( m_rand_uint32 >> 13 );
}

int np02_shape_test::execute(){
int error_code = 0;
const time_t start_time = time(NULL);
m_rand_uint32 = m_shape_test_rand_seed ^ 4294967291ul;

switch( m_shape_test_number){
    case 1: error_code = execute_test_1(); break;
    case 2: error_code = execute_test_2(); break;
    case 3: error_code = execute_test_3(); break;
    case 4: error_code = execute_test_4(); break;
    case 5: error_code = execute_test_5(); break;
    case 0:
    default:
        std::cout << "attempt to run shape test " << m_shape_test_number 
            << " fails.  no such test\n";
        break;
    }

const time_t done_time = time(NULL);
std::cout << "done  " << ctime(&done_time);
std::cout << "elapsed_time=" << (done_time-start_time) << " s\n";
return error_code;
}

void np02_shape_test::make_rand_shapes( np02_shp_alloc *shp_alloc,
    const int& ww, const int& hh,  const double& basic_length,
    std::vector<np02_shape *> *shapes,
    std::vector<np02_bmp_color> *colors ){
for(int i = 0; i < ww; ++i){
    const double x_ctr = (static_cast<double>(i)+0.5) * basic_length;
    for(int j = 0; j < hh; ++j){
        /* create random shape & color */
        const double y_ctr = (static_cast<double>(j)+0.5) * basic_length;
        np02_bmp_color color(0,0,0,false);
        np02_shape *shape = make_rand_shape( shp_alloc,
            np02_xy(x_ctr, y_ctr), basic_length, &color );
        shapes->push_back(shape);
        colors->push_back(color);
        }
    }
}


np02_shape *np02_shape_test::make_rand_shape( np02_shp_alloc *shp_alloc,
    const np02_xy& shp_ctr, const double& basic_length,
    np02_bmp_color *color ){
np02_shape *shape = NULL;
advance_rand();
double x0 = get_rand_dbl(-basic_length*0.75, basic_length*0.75);
advance_rand();
double y0 = get_rand_dbl(-basic_length*0.75, basic_length*0.75);
advance_rand();
double x1 = get_rand_dbl(-basic_length*0.75, basic_length*0.75);
advance_rand();
double y1 = get_rand_dbl(-basic_length*0.75, basic_length*0.75);
advance_rand();
double r = get_rand_dbl(0.0, basic_length/2.0);
advance_rand();
double rot_deg = static_cast<double>(get_rand_int()%24)*15.0;
np02_bmp_color clr(0,0,0);
switch(get_rand_int()%9){
    default:
    case 0:
    case 3:
        {
        np02_circle *circle = (NULL == shp_alloc) ? 
            new np02_circle() : shp_alloc->alloc_circle();
        circle->init(np02_xy(shp_ctr.get_x()+x0,shp_ctr.get_y()+y0),r);
        shape = circle;

        /* create random color */
        advance_rand();        
        switch(get_rand_int()%4){
            default:
            case  0: clr=np02_bmp_color::bl_olive  ();break;
            case  1: clr=np02_bmp_color::bl_lime   ();break;
            case  2: clr=np02_bmp_color::bl_green  ();break;
            }
        }
        break;
    case 1:
    case 4:
    case 6:
    case 7:
        {
        np02_rect *rect = (NULL == shp_alloc) ? 
            new np02_rect() : shp_alloc->alloc_rect();
        rect->init(np02_xy(shp_ctr.get_x()+x0, shp_ctr.get_y()+y0),
            fabs(x1), fabs(y1), rot_deg);
        shape = rect;

        /* create random color */
        advance_rand();        
        switch(get_rand_int()%2){
            default:
            case  0: clr=np02_bmp_color::bl_blue   ();break;
            case  1: clr=np02_bmp_color::bl_navy   ();break;
            }
        }
        break;
    case 2:
    case 5:
    case 8:
        {
        np02_line_seg *line_seg = (NULL == shp_alloc) ? 
            new np02_line_seg() : shp_alloc->alloc_line_seg();
        line_seg->init(np02_xy(shp_ctr.get_x()+x0,shp_ctr.get_y()+y0), 
            np02_xy(shp_ctr.get_x()+x1,shp_ctr.get_y()+y1), r);
        shape = line_seg;

        /* create random color */
        advance_rand();

        switch(get_rand_int()%3){
            default:
            case 0: clr=np02_bmp_color::bl_aqua ();break;
            case 1: clr=np02_bmp_color::bl_teal();break;
            case 2: clr=np02_bmp_color::bl_purple ();break;
            }
        }
        break;
    }

if(NULL != color){ *color = clr; }
return shape;
}

void np02_shape_test::free_shapes( np02_shp_alloc *shp_alloc,
    std::vector<np02_shape *> *shapes ){
int e = 0;
if(NULL != shp_alloc){
    e = shp_alloc->verify_data(AA_ERR_BUF(),
        AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR());
    AA_ALWAYS_ASSERT(0 == e);
    }

if (NULL != shapes){
    std::vector<np02_shape *>::const_iterator shp_itr = shapes->begin();
    for(; shp_itr != shapes->end(); ++shp_itr){
        np02_shape *shape = *shp_itr;
        if(NULL == shp_alloc){
            delete shape;
            }
        else{
            shp_alloc->free_shape(shape);
            }
        }
    shapes->clear();
    }

if(NULL != shp_alloc){
    e = shp_alloc->verify_data(AA_ERR_BUF(),
        AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR());
    AA_ALWAYS_ASSERT(0 == e);
    }
}

int np02_shape_test::execute_test_1(){
AA_INCR_CALL_DEPTH();
int err_cnt = 0;
std::cout << "shape test 1\n";
std::cout << "iterations: " << m_shape_test_iteration_count;
std::cout << "  rand_seed: " << m_shape_test_rand_seed << "\n";

m_iteration = 0;
for(; m_iteration < m_shape_test_iteration_count; ++m_iteration ){
    switch(m_iteration){
        case 0: err_cnt += execute_test_1_0(); break;
        case 1: err_cnt += execute_test_1_1(); break;
        case 2: err_cnt += execute_test_1_2(); break;
        default: err_cnt += execute_test_1_i(); break;
        }
    }
std::cout << "errors: " << err_cnt << "\n";
if( err_cnt > 0 ){ std::cout << AA_ERR_BUF(); }
AA_DECR_CALL_DEPTH();
return err_cnt;
}

int np02_shape_test::execute_test_1_0(){
AA_INCR_CALL_DEPTH();
int err_cnt = 0;
int e = 0;

std::vector<np02_shape *> shape_vec;
std::vector<np02_bmp_color> color_vec;


np02_circle circle;
circle.init(np02_xy(10.0, 10.0), 7.0);
shape_vec.push_back(&circle);
color_vec.push_back(np02_bmp_color(255, 128, 0, false));
e = circle.verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR());
AA_ALWAYS_ASSERT(0 == e);
err_cnt += e;

np02_line_seg line_seg;
line_seg.init(np02_xy(-8.0, -11.0), np02_xy(-5.0, 8.0), 5.0);
shape_vec.push_back(&line_seg);
color_vec.push_back(np02_bmp_color(128, 0, 255, true));
e = line_seg.verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR());
AA_ALWAYS_ASSERT(0 == e);
err_cnt += e;


np02_line_seg line_seg2;
line_seg2.init(np02_xy(-16.0, -2.0), np02_xy(-3.0, 4.0), 3.0);
shape_vec.push_back(&line_seg2);
color_vec.push_back(np02_bmp_color(128, 0, 255, true));
e = line_seg2.verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR());
AA_ALWAYS_ASSERT(0 == e);
err_cnt += e;


np02_line_seg line_seg3;
line_seg3.init(np02_xy(-17.0, 12.0), np02_xy(-14.0, 13.0), 2.0);
shape_vec.push_back(&line_seg3);
color_vec.push_back(np02_bmp_color(128, 0, 255, true));
e = line_seg3.verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR());
AA_ALWAYS_ASSERT(0 == e);
err_cnt += e;

np02_rect rect;
rect.init(np02_xy(8.0, -10.0), 11.0, 7.0, 30.0);
shape_vec.push_back(&rect);
color_vec.push_back(np02_bmp_color(0, 255, 128, false));
e = rect.verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR());
AA_ALWAYS_ASSERT(0 == e);
err_cnt += e;

np02_circle circle2;
circle2.init(np02_xy(17.0, 1.0), 2.0);
shape_vec.push_back(&circle2);
color_vec.push_back(np02_bmp_color(255, 128, 0, false));
e = circle2.verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR());
AA_ALWAYS_ASSERT(0 == e);
err_cnt += e;

np02_rect rect2;
rect2.init(np02_xy(18.0, -11.0), 2.0, 1.5, -10.0);
shape_vec.push_back(&rect2);
color_vec.push_back(np02_bmp_color(0, 255, 128, false));
e = rect2.verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR());
AA_ALWAYS_ASSERT(0 == e);
err_cnt += e;


np02_rect rect3;
rect3.init(np02_xy(5.0, -19.0), 3.0, 2.0, -10.0);
shape_vec.push_back(&rect3);
color_vec.push_back(np02_bmp_color(0, 255, 128, false));
e = rect3.verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR());
AA_ALWAYS_ASSERT(0 == e);
err_cnt += e;

np02_circle circle5;
circle5.init(np02_xy(-1.0, -9.0), 2.0);
shape_vec.push_back(&circle5);
color_vec.push_back(np02_bmp_color(255, 128, 0, false));
e = circle5.verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR());
AA_ALWAYS_ASSERT(0 == e);
err_cnt += e;


np02_bmp_file_init_params bmp_init_params;
const double pixel_num = 8.0;
bmp_init_params.width_px = 320;
bmp_init_params.height_px = 320;
np02_bmp_file bmp_file(bmp_init_params);
const np02_xy xy_min(-20.0,-20.0);
size_t i,j;
for( i = 0; i < shape_vec.size(); ++i){
    shape_vec.at(i)->write_bmp_file(xy_min, pixel_num, color_vec.at(i), &bmp_file);
    }
for( i = 0; i < shape_vec.size(); ++i){
    const np02_shape *shape_i = shape_vec.at(i);
    for( j = i+1; j < shape_vec.size(); ++j){
        const np02_shape *shape_j = shape_vec.at(j);
        np02_xy xy_a, xy_b;
        const double d = shape_i->get_distance_from_shape(shape_j,
            &xy_a, &xy_b );
        const int32_t a_ii = static_cast<int32_t>(
            (xy_a.get_x() - xy_min.get_x()) * pixel_num);
        const int32_t a_jj = static_cast<int32_t>(
            (xy_a.get_y() - xy_min.get_y()) * pixel_num);
        const int32_t b_ii = static_cast<int32_t>(
            (xy_b.get_x() - xy_min.get_x()) * pixel_num);
        const int32_t b_jj = static_cast<int32_t>(
            (xy_b.get_y() - xy_min.get_y()) * pixel_num);
        if( d < 0.0 ){
            np02_xy xy_ab_ave( (xy_a.get_x() + xy_b.get_x())/2.0,
                (xy_a.get_y() + xy_b.get_y())/2.0 );
            const double ii_ctr =
                ( xy_ab_ave.get_x() - xy_min.get_x() ) * pixel_num;
            const double jj_ctr =
                ( xy_ab_ave.get_y() - xy_min.get_y() ) * pixel_num;
            const double rr = fabs(d) * pixel_num;

            bmp_file.draw_circle(ii_ctr,jj_ctr,rr,np02_bmp_color(255,0,0, true));
            }
        bmp_file.draw_line(a_ii, a_jj, b_ii, b_jj,
            np02_bmp_color(255,255,255,true));
        }
    }
bmp_file.write_file( "shape_test_0.bmp" );

AA_DECR_CALL_DEPTH();
return err_cnt;
}

int np02_shape_test::execute_test_1_1(){
int err_cnt = 0;
int e = 0;
AA_INCR_CALL_DEPTH();
err_cnt += execute_test_1_i();
AA_DECR_CALL_DEPTH();
return err_cnt;
}

int np02_shape_test::execute_test_1_2(){
int err_cnt = 0;
AA_INCR_CALL_DEPTH();
err_cnt += execute_test_1_i();
AA_DECR_CALL_DEPTH();
return err_cnt;
}

int np02_shape_test::execute_test_1_i(){
AA_INCR_CALL_DEPTH();
int e = 0;
int err_cnt = 0;
size_t i,j;
std::vector<np02_shape *> shape_vec;
std::vector<np02_bmp_color> color_vec;
np02_shp_alloc *shp_alloc = NULL;

advance_rand();
int ww = 3 + (get_rand_int() % 6);
advance_rand();
int hh = 3 + (get_rand_int() % 6);
advance_rand();
double basic_length = 1.0;
switch(get_rand_int() % 7){
    default: 
    case 0:  basic_length = 0.01; break;
    case 1:  basic_length = 0.1; break;
    case 2:  basic_length = 1.0; break;
    case 3:  basic_length = 10.0; break;
    case 4:  basic_length = 100.0; break;
    case 5:  basic_length = 1000.0; break;
    case 6:  basic_length = 10000.0; break;
    }
advance_rand();
if(50 > (get_rand_int() % 100)){
    shp_alloc = new np02_shp_alloc();
    e = shp_alloc->verify_data(AA_ERR_BUF(),
        AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR());
    AA_ALWAYS_ASSERT(0 == e);
    }
make_rand_shapes( shp_alloc, ww, hh, basic_length, &shape_vec, &color_vec );

bool should_write_file = false;
if( ( m_iteration < 10 ) ||
    ( ( m_iteration < 100 ) && ( (m_iteration % 10) == 0 ) ) ||
    ( ( m_iteration < 1000 ) && ( (m_iteration % 100) == 0 ) ) ||
    ( ( m_iteration < 10000 ) && ( (m_iteration % 1000) == 0 ) ) ){
    should_write_file = true;
    }

const np02_xy xy_min(0.0,0.0);

bool err_xy_found = false;
np02_xy err_xy_a(0.0, 0.0);
np02_xy err_xy_b(0.0, 0.0);

for( i = 0; i < shape_vec.size(); ++i){
    e = shape_vec.at(i)->verify_data(AA_ERR_BUF(),
        AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR());
    AA_ALWAYS_ASSERT(0 == e);
    err_cnt += e;
    if((e > 0) && !err_xy_found ){
        err_xy_found = true;
        should_write_file = true;
        shape_vec.at(i)->get_bb(&err_xy_a, &err_xy_b );
        }
    }

for( i = 0; i < shape_vec.size(); ++i){
    const np02_shape *shape_i = shape_vec.at(i);
    for( j = i+1; j < shape_vec.size(); ++j){
        const np02_shape *shape_j = shape_vec.at(j);
        np02_xy xy_a, xy_b;
        const double d = shape_i->get_distance_from_shape(shape_j,
            &xy_a, &xy_b );
        e = np02_shape::verify_distance_from_shape_result(
            shape_i, shape_j, xy_a, xy_b, d, AA_ERR_BUF(),
            AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR());
        AA_ALWAYS_ASSERT(0 == e);
        err_cnt += e;
        if((e > 0) && !err_xy_found ){
            err_xy_found = true;
            should_write_file = true;
            err_xy_a=xy_a;
            err_xy_b=xy_b;
            }
        }
    }

if(should_write_file){
    np02_bmp_file_init_params bmp_init_params;
    static const double pixels_per_basic_length = 256.0;
    const double pixel_num = pixels_per_basic_length/basic_length;
    bmp_init_params.width_px = ww * pixels_per_basic_length;
    bmp_init_params.height_px = hh * pixels_per_basic_length;
    np02_bmp_file bmp_file(bmp_init_params);
    char bmp_file_name[255];
    sprintf(&(bmp_file_name[0]), "shape_test1_%i.bmp", m_iteration );

    for( i = 0; i < shape_vec.size(); ++i){
        shape_vec.at(i)->write_bmp_file(xy_min, pixel_num,
            color_vec.at(i), &bmp_file);
        }

    for( i = 0; i < shape_vec.size(); ++i){
        const np02_shape *shape_i = shape_vec.at(i);
        for( j = i+1; j < shape_vec.size(); ++j){
            const np02_shape *shape_j = shape_vec.at(j);
            np02_xy xy_a, xy_b;
    
            const double d = shape_i->get_distance_from_shape(shape_j,
                &xy_a, &xy_b );
    
            const int32_t a_ii = static_cast<int32_t>(
                (xy_a.get_x() - xy_min.get_x()) * pixel_num);
            const int32_t a_jj = static_cast<int32_t>(
                (xy_a.get_y() - xy_min.get_y()) * pixel_num);
            const int32_t b_ii = static_cast<int32_t>(
                (xy_b.get_x() - xy_min.get_x()) * pixel_num);
            const int32_t b_jj = static_cast<int32_t>(
                (xy_b.get_y() - xy_min.get_y()) * pixel_num);
            if(fabs(d) < (1.25 * basic_length) ){
                if( d < 0.0 ){
                    np02_xy xy_ab_ave( (xy_a.get_x() + xy_b.get_x())/2.0,
                        (xy_a.get_y() + xy_b.get_y())/2.0 );
                    const double ii_ctr =
                        ( xy_ab_ave.get_x() - xy_min.get_x() ) * pixel_num;
                    const double jj_ctr =
                        ( xy_ab_ave.get_y() - xy_min.get_y() ) * pixel_num;
                    const double rr = fabs(d) * pixel_num;
        
                    bmp_file.draw_circle(ii_ctr,jj_ctr,rr,np02_bmp_color(255,0,0, true));
                    }
                bmp_file.draw_line(a_ii, a_jj, b_ii, b_jj,
                    np02_bmp_color(255,255,255,true));
                }
            }
        }

    if(err_xy_found){
        const np02_xy text_pos((err_xy_a.get_x() + err_xy_b.get_x())/2.0,
                                   (err_xy_a.get_y() + err_xy_b.get_y())/2.0);
        const double d_err_ab = err_xy_a.get_distance_to(err_xy_b);
        double hh = d_err_ab * pixel_num;
        if(hh < 10){ hh = 10.0; }
        const double ii_text =
            ( text_pos.get_x() - xy_min.get_x() ) * pixel_num;
        const double jj_text =
            ( text_pos.get_y() - xy_min.get_y() ) * pixel_num;
        std::cout << __FILE__ << "[" << __LINE__ << "] "
            << "shape error  itr=" << m_iteration
            << "  drawing=" << bmp_file_name 
            << "  xy(" << text_pos.get_x() << "," << text_pos.get_y() << ")"
            << "  pixel(" << ii_text << "," << jj_text << ")\n";
        bmp_file.draw_text("Error", ii_text, jj_text, hh, 15.0, np02_bmp_color::white());
        }

    std::string prog_str = np02_test_main::get_prog_name() + std::string(" ") + np02_test_main::get_version_str();
    bmp_file.draw_text(prog_str.c_str(), 25.0, bmp_init_params.height_px-75.0,
        50.0, 0.0, np02_bmp_color::gray());
    bmp_file.draw_text(&(bmp_file_name[0]), 25.0, 25.0, 50.0, 0.0,
        np02_bmp_color::gray());

    bmp_file.write_file( &(bmp_file_name[0]) );
    }

free_shapes( shp_alloc, &shape_vec );
if( NULL != shp_alloc ){ delete shp_alloc; }

AA_DECR_CALL_DEPTH();

return err_cnt;
}


int np02_shape_test::execute_test_2(){ 
AA_INCR_CALL_DEPTH();
int err_cnt = 0;
std::cout << "shape test 2\n";
std::cout << "iterations: " << m_shape_test_iteration_count;
std::cout << "  rand_seed: " << m_shape_test_rand_seed << "\n";
m_iteration = 0;
for(; m_iteration < m_shape_test_iteration_count; ++m_iteration ){
    switch(m_iteration % 2){
        default:
        case 0: err_cnt += execute_test_2_0(); break;
        case 1: err_cnt += execute_test_2_1(); break;
        }
    }
std::cout << "errors: " << err_cnt << "\n";
if( err_cnt > 0 ){ std::cout << AA_ERR_BUF(); }
AA_DECR_CALL_DEPTH();
return err_cnt;
}


int np02_shape_test::execute_test_2_0(){ 
AA_INCR_CALL_DEPTH();
int e = 0;
int err_cnt = 0;
int local_err_cnt = 0;
size_t i,j;
std::vector<np02_shape *> shape_vec;
std::vector<np02_bmp_color> color_vec;
np02_shp_alloc *shp_alloc = NULL;

/* make random shapes */
advance_rand();
int ww = 3 + (get_rand_int() % 6);
advance_rand();
int hh = 3 + (get_rand_int() % 6);
advance_rand();
double basic_length = 1.0;
switch(get_rand_int() % 7){
    default: 
    case 0:  basic_length = 0.01; break;
    case 1:  basic_length = 0.1; break;
    case 2:  basic_length = 1.0; break;
    case 3:  basic_length = 10.0; break;
    case 4:  basic_length = 100.0; break;
    case 5:  basic_length = 1000.0; break;
    case 6:  basic_length = 10000.0; break;
    }
if(50 > (get_rand_int() % 100)){
    shp_alloc = new np02_shp_alloc();
    e = shp_alloc->verify_data(AA_ERR_BUF(),
        AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR());
    AA_ALWAYS_ASSERT(0 == e);
    }
make_rand_shapes( shp_alloc, ww, hh, basic_length, &shape_vec, &color_vec );
const np02_xy xy_min(0.0,0.0);

bool should_write_file = false;
if( ( m_iteration < 10 ) ||
    ( ( m_iteration < 100 ) && ( (m_iteration % 10) == 0 ) ) ||
    ( ( m_iteration < 1000 ) && ( (m_iteration % 100) == 0 ) ) ||
    ( ( m_iteration < 10000 ) && ( (m_iteration % 1000) == 0 ) ) ){
    should_write_file = true;
    }

/* make randomized locator grid */
advance_rand();
double extra_search_d = 0.0;
switch(get_rand_int() % 8){
    case 0:
        extra_search_d = 0.0;
        break;
    case 1:
        advance_rand();
        extra_search_d = get_rand_dbl( basic_length, 2.0*basic_length );
        break;
    default:
        advance_rand();
        extra_search_d = get_rand_dbl( 0.0, basic_length );
        break;
    }
const double sq_sz = get_rand_dbl(basic_length/4.0, basic_length);
np02_loc_grid_dim loc_grid_dim;
uint16_t w = static_cast<uint16_t>((ww*basic_length)/sq_sz);
if( w > 3 ){ advance_rand();  w -= (get_rand_int() % 3); }
uint16_t h = static_cast<uint16_t>((hh*basic_length)/sq_sz);
if( h > 3 ){ advance_rand();  h -= (get_rand_int() % 3); }
loc_grid_dim.set_w(w);
loc_grid_dim.set_h(h);
advance_rand();
loc_grid_dim.set_x_min(get_rand_dbl(0.0, basic_length));
advance_rand();
loc_grid_dim.set_y_min(get_rand_dbl(0.0, basic_length));
loc_grid_dim.set_sq_size(sq_sz);

np02_loc_grid loc_grid_local;
np02_loc_grid *loc_grid;

if( NULL != shp_alloc ){
    loc_grid = shp_alloc->alloc_loc_grid();
    //std::cout << "loc_grid = shp_alloc->alloc_loc_grid()" << "\n";
    }
else if(50 > (get_rand_int() % 100)){
    loc_grid = new np02_loc_grid();
    //std::cout << "loc_grid = new np02_loc_grid()" << "\n";
    }
else{
    loc_grid = &loc_grid_local;
    //std::cout << "loc_grid = &loc_grid_local" << "\n";
    }

e = loc_grid->verify_data(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR());
AA_ALWAYS_ASSERT(0 == e);

loc_grid->set_extra_search_d(extra_search_d);
loc_grid->init_loc_grid(loc_grid_dim);

e = loc_grid->verify_data(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR());
AA_ALWAYS_ASSERT(0 == e);

/* Put shapes in locator grid */
bool err_xy_found = false;
np02_xy err_xy_a(0.0, 0.0);
np02_xy err_xy_b(0.0, 0.0);
np02_shape *shape = NULL;
for( i = 0; i < shape_vec.size(); ++i){
    shape = shape_vec.at(i);

    e = shape->verify_data(AA_ERR_BUF(),
        AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR());
    AA_ALWAYS_ASSERT(0 == e);
    err_cnt += e;
    if((e > 0) && !err_xy_found ){
        err_xy_found = true;
        should_write_file = true;
        shape->get_bb(&err_xy_a, &err_xy_b );
        }

    loc_grid->insert_shape_in_loc_grid(shape);
    
    e = shape->verify_data(AA_ERR_BUF(),
        AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR());
    AA_ALWAYS_ASSERT(0 == e);
    err_cnt += e;
    if((e > 0) && !err_xy_found ){
        err_xy_found = true;
        should_write_file = true;
        shape->get_bb(&err_xy_a, &err_xy_b );
        }
    }

e = loc_grid->verify_data(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR());
AA_ALWAYS_ASSERT(0 == e);
err_cnt += e;

/* 
From each shape,
  find nearest distance to every other shape
  find all local shapes with locator grid
  check that all shapes within locator grid extra-search-distance are found
    in the local locator grid search.
*/
np02_loc_grid::shape_vec local_shapes;
for( i = 0; i < shape_vec.size(); ++i){
    const np02_shape *shape_i = shape_vec.at(i);
    loc_grid->get_shapes_near_shape(shape_i, &local_shapes);
    std::sort(local_shapes.begin(), local_shapes.end());
    for( j = 0; j < shape_vec.size(); ++j){
        if( i != j ){
            const np02_shape *shape_j = shape_vec.at(j);
            np02_xy xy_a, xy_b;
            const double d = shape_i->get_distance_from_shape(shape_j,
                &xy_a, &xy_b );
            if( d <= extra_search_d ){
                np02_loc_grid::shape_vec_citr search_itr =
                    std::lower_bound( local_shapes.begin(), local_shapes.end(),
                        const_cast<np02_shape *>(shape_j));
                if( (search_itr == local_shapes.end()) || 
                    (shape_j != *search_itr)){
                    AA_ALWAYS_ASSERT(false);
                    if(!err_xy_found ){
                        err_xy_found = true;
                        should_write_file = true;
                        shape_i->get_bb(&err_xy_a, &err_xy_b );
                        }
                    }
                }
            }
        }
    local_shapes.clear();
    }

if(should_write_file){
    np02_bmp_file_init_params bmp_init_params;
    static const double pixels_per_basic_length = 256.0;
    const double pixel_num = pixels_per_basic_length/basic_length;
    bmp_init_params.width_px = static_cast<int32_t>(ww*pixels_per_basic_length);
    bmp_init_params.height_px =
        static_cast<int32_t>(hh * pixels_per_basic_length);
    np02_bmp_file bmp_file(bmp_init_params);
    char bmp_file_name[255];
    sprintf(&(bmp_file_name[0]), "shape_test2_%i.bmp", m_iteration );

    for( i = 0; i < shape_vec.size(); ++i){
        shape_vec.at(i)->write_bmp_file(xy_min, pixel_num,
            color_vec.at(i), &bmp_file);
        }

    loc_grid->write_bmp_file(xy_min, pixel_num, np02_bmp_color::gray(),
        &bmp_file);

    if(err_xy_found){
        const np02_xy text_pos((err_xy_a.get_x() + err_xy_b.get_x())/2.0,
                                   (err_xy_a.get_y() + err_xy_b.get_y())/2.0);
        const double d_err_ab = err_xy_a.get_distance_to(err_xy_b);
        double hh = d_err_ab * pixel_num;
        if(hh < 10){ hh = 10.0; }
        const double ii_text =
            ( text_pos.get_x() - xy_min.get_x() ) * pixel_num;
        const double jj_text =
            ( text_pos.get_y() - xy_min.get_y() ) * pixel_num;
        std::cout << __FILE__ << "[" << __LINE__ << "] "
            << "shape error  itr=" << m_iteration
            << "  drawing=" << bmp_file_name 
            << "  xy(" << text_pos.get_x() << "," << text_pos.get_y() << ")"
            << "  pixel(" << ii_text << "," << jj_text << ")\n";
        bmp_file.draw_text("Error", ii_text, jj_text, hh, 15.0, np02_bmp_color::white());
        }

    std::string prog_str = np02_test_main::get_prog_name() + std::string(" ") +
        np02_test_main::get_version_str();
    bmp_file.draw_text(prog_str.c_str(), 25.0, bmp_init_params.height_px-75.0,
        50.0, 0.0, np02_bmp_color::gray());
    bmp_file.draw_text(&(bmp_file_name[0]), 25.0, 25.0, 50.0, 0.0,
        np02_bmp_color::gray());

    bmp_file.write_file( &(bmp_file_name[0]) );
    }

err_xy_found = false;

/* 
Move some shapes
*/
for( i = 0; i < shape_vec.size(); ++i){
    np02_shape *shape_i = shape_vec.at(i);

    local_err_cnt = 0;
    e = shape_i->verify_data(AA_ERR_BUF(),
        AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR());
    AA_ALWAYS_ASSERT(0 == e);
    local_err_cnt += e;
    e = (shape_i->get_loc_grid() == loc_grid)? 0 : 1;
    AA_ALWAYS_ASSERT(0 == e);
    local_err_cnt += e;

    np02_xy shape_i_bb_xy_min, shape_i_bb_xy_max;
    shape_i->get_bb(&shape_i_bb_xy_min, &shape_i_bb_xy_max);
    advance_rand();
    const double rot_ctr_x = get_rand_dbl(shape_i_bb_xy_min.get_x(),
        shape_i_bb_xy_max.get_x());
    advance_rand();
    const double rot_ctr_y = get_rand_dbl(shape_i_bb_xy_min.get_y(),
        shape_i_bb_xy_max.get_y());
    advance_rand();
    const double rot_deg = static_cast<double>(get_rand_int() % 721)-360.0;
    const double translation_dx = get_rand_dbl(-basic_length, basic_length);
    advance_rand();
    const double translation_dy = get_rand_dbl(-basic_length, basic_length);
    advance_rand();
    switch (get_rand_int() % 10) {
        default:
        case 0:
        case 1:
            /* don't move */
            break;
        case 2:
        case 3:
            /* rotate, automatically update locator grid */
            shape->rotate(np02_xy(rot_ctr_x, rot_ctr_y), rot_deg);
            break;
        case 4:
        case 5:
            /* translate, automatically update locator grid */
            shape->translate(np02_xy(translation_dx, translation_dy));
            break;
        case 6:
        case 7:
        case 8:
            {
            /* rotate, explicitly update locator grid */
            loc_grid->remove_shape_from_loc_grid(shape_i);

            e = shape_i->verify_data(AA_ERR_BUF(),
                AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR());
            AA_ALWAYS_ASSERT(0 == e);
            local_err_cnt += e;
            e = (shape_i->get_loc_grid() == NULL)? 0 : 1;
            AA_ALWAYS_ASSERT(0 == e);
            local_err_cnt += e;

            shape_i->rotate_no_loc_grid(
                np02_xy(rot_ctr_x, rot_ctr_y), rot_deg);

            e = shape_i->verify_data(AA_ERR_BUF(),
                AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR());
            AA_ALWAYS_ASSERT(0 == e);
            local_err_cnt += e;

            loc_grid->insert_shape_in_loc_grid(shape_i);
            }
            break;
        case 9:
        case 10:
        case 11:
            {
            /* translate, explicitly update locator grid */
            loc_grid->remove_shape_from_loc_grid(shape_i);

            e = shape_i->verify_data(AA_ERR_BUF(),
                AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR());
            AA_ALWAYS_ASSERT(0 == e);
            local_err_cnt += e;
            e = (shape_i->get_loc_grid() == NULL)? 0 : 1;
            AA_ALWAYS_ASSERT(0 == e);
            local_err_cnt += e;

            shape_i->translate_no_loc_grid(
                np02_xy(translation_dx, translation_dy));

            e = shape_i->verify_data(AA_ERR_BUF(),
                AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR());
            AA_ALWAYS_ASSERT(0 == e);
            local_err_cnt += e;

            loc_grid->insert_shape_in_loc_grid(shape_i);
            }
            break;
        }

    e = shape_i->verify_data(AA_ERR_BUF(),
        AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR());
    AA_ALWAYS_ASSERT(0 == e);
    local_err_cnt += e;
    e = (shape_i->get_loc_grid() == loc_grid)? 0 : 1;
    AA_ALWAYS_ASSERT(0 == e);
    local_err_cnt += e;

    err_cnt += local_err_cnt;
    if((local_err_cnt > 0) && !err_xy_found ){
        err_xy_found = true;
        should_write_file = true;
        shape_i->get_bb(&err_xy_a, &err_xy_b );
        }
    }

/* 
From each shape,
  find nearest distance to every other shape
  find all local shapes with locator grid
  check that all shapes within locator grid extra-search-distance are found
    in the local locator grid search.
*/
for( i = 0; i < shape_vec.size(); ++i){
    const np02_shape *shape_i = shape_vec.at(i);
    loc_grid->get_shapes_near_shape(shape_i, &local_shapes);
    std::sort(local_shapes.begin(), local_shapes.end());
    for( j = 0; j < shape_vec.size(); ++j){
        if( i != j ){
            const np02_shape *shape_j = shape_vec.at(j);
            np02_xy xy_a, xy_b;
            const double d = shape_i->get_distance_from_shape(shape_j,
                &xy_a, &xy_b );
            if( d <= extra_search_d ){
                np02_loc_grid::shape_vec_citr search_itr =
                    std::lower_bound( local_shapes.begin(), local_shapes.end(),
                        const_cast<np02_shape *>(shape_j));
                if( (search_itr == local_shapes.end()) || 
                    (shape_j != *search_itr)){
                    AA_ALWAYS_ASSERT(false);
                    if(!err_xy_found ){
                        err_xy_found = true;
                        should_write_file = true;
                        shape_i->get_bb(&err_xy_a, &err_xy_b );
                        }
                    }
                }
            }
        }
    local_shapes.clear();
    }

if(should_write_file){
    np02_bmp_file_init_params bmp_init_params;
    static const double pixels_per_basic_length = 256.0;
    const double pixel_num = pixels_per_basic_length/basic_length;
    bmp_init_params.width_px =
        static_cast<int32_t>(ww * pixels_per_basic_length);
    bmp_init_params.height_px =
        static_cast<int32_t>(hh * pixels_per_basic_length);
    np02_bmp_file bmp_file(bmp_init_params);
    char bmp_file_name[255];
    sprintf(&(bmp_file_name[0]), "shape_test2a_%i.bmp", m_iteration );

    for( i = 0; i < shape_vec.size(); ++i){
        shape_vec.at(i)->write_bmp_file(xy_min, pixel_num,
            color_vec.at(i), &bmp_file);
        }

    loc_grid->write_bmp_file(xy_min, pixel_num, np02_bmp_color::gray(),
        &bmp_file);

    if(err_xy_found){
        const np02_xy text_pos((err_xy_a.get_x() + err_xy_b.get_x())/2.0,
                                   (err_xy_a.get_y() + err_xy_b.get_y())/2.0);
        const double d_err_ab = err_xy_a.get_distance_to(err_xy_b);
        double hh = d_err_ab * pixel_num;
        if(hh < 10){ hh = 10.0; }
        const double ii_text =
            ( text_pos.get_x() - xy_min.get_x() ) * pixel_num;
        const double jj_text =
            ( text_pos.get_y() - xy_min.get_y() ) * pixel_num;
        std::cout << __FILE__ << "[" << __LINE__ << "] "
            << "shape error  itr=" << m_iteration
            << "  drawing=" << bmp_file_name 
            << "  xy(" << text_pos.get_x() << "," << text_pos.get_y() << ")"
            << "  pixel(" << ii_text << "," << jj_text << ")\n";
        bmp_file.draw_text("Error", ii_text, jj_text, hh, 15.0, np02_bmp_color::white());
        }

    std::string prog_str = np02_test_main::get_prog_name() + std::string(" ") +
        np02_test_main::get_version_str();
    bmp_file.draw_text(prog_str.c_str(), 25.0, bmp_init_params.height_px-75.0,
        50.0, 0.0, np02_bmp_color::gray());
    bmp_file.draw_text(&(bmp_file_name[0]), 25.0, 25.0, 50.0, 0.0,
        np02_bmp_color::gray());

    bmp_file.write_file( &(bmp_file_name[0]) );
    }

e = loc_grid->verify_data(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR());
AA_ALWAYS_ASSERT(0 == e);
err_cnt += e;

free_shapes( shp_alloc, &shape_vec );

e = loc_grid->verify_data(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR());
AA_ALWAYS_ASSERT(0 == e);
err_cnt += e;

if( NULL != shp_alloc ){ delete shp_alloc; }
loc_grid = NULL;

AA_DECR_CALL_DEPTH();

return err_cnt;
}

int np02_shape_test::execute_test_2_1(){ return execute_test_2_0(); }

int np02_shape_test::execute_test_3(){ return 0; }
int np02_shape_test::execute_test_4(){ return 0; }
int np02_shape_test::execute_test_5(){ return 0; }



} /* namespace np02 */
