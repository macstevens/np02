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
#include <sstream>
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

const double np02_shape::m_pi = 3.1415926535897932384626433832795;
const double np02_shape::m_degrees_per_radian
    = 180.0 / 3.1415926535897932384626433832795;
const double np02_shape::m_little_ratio
    = 0.0000152587890625; /* 0.5 ^ 16 */
const double np02_shape::m_little_ratio_sq =
    ( np02_shape::m_little_ratio * np02_shape::m_little_ratio );
const double np02_shape::m_small_ratio 
    = 0.00000000023283064365386962890625; /* 0.5 ^ 32 */
const double np02_shape::m_small_ratio_sq =
    ( np02_shape::m_small_ratio * np02_shape::m_small_ratio );
const double np02_shape::m_tiny_ratio
    = 5.4210108624275221700372640043497e-20; /* 0.5 ^ 64 */
const double np02_shape::m_teensy_ratio
    = 2.9387358770557187699218413430556e-39; /* 0.5 ^ 128 */
const double np02_shape::m_teensy_ratio_sq =
    ( np02_shape::m_teensy_ratio * np02_shape::m_teensy_ratio );
const double np02_shape::m_googol = 1e100;


void np02_xy::rotate(const np02_xy& rot_ctr, const double& rot_deg){
    const np02_xy rot_arm_initial(m_x - rot_ctr.get_x(),
        m_y - rot_ctr.get_y());
    if( (rot_arm_initial.get_x() != 0.0) || 
        (rot_arm_initial.get_y() != 0.0) ){
        double cos_rot = 1.0;
        double sin_rot = 0.0;
        np02_shape::cos_sin_rot_deg( rot_deg, &cos_rot, &sin_rot );
        const np02_xy rot_arm_final(
            (rot_arm_initial.get_x() * cos_rot) -
            (rot_arm_initial.get_y() * sin_rot),
            (rot_arm_initial.get_x() * sin_rot) +
            (rot_arm_initial.get_y() * cos_rot));
        set_x( rot_ctr.get_x() + rot_arm_final.get_x() );
        set_y( rot_ctr.get_y() + rot_arm_final.get_y() );
        }
    }

uint64_t np02_xy::hash( const uint64_t& h_in ) const{
uint64_t h = h_in;
#if defined( CF01_SUPPORT )
h = cf01_obj_hash( h, m_x );
h = cf01_obj_hash( h, m_y );
#else
h += reinterpret_cast<cf01_uint64>(m_x);
h ^= ((h << 37) | (h >> 27));
h += reinterpret_cast<cf01_uint64>(m_y);
h ^= ((h << 47) | (h >> 17));
#endif
return h;
}


np02_shape_owner::np02_shape_owner(){
}

np02_shape_owner::~np02_shape_owner(){
/* derived class should delete shape */
}

uint64_t np02_shape_owner::hash( const uint64_t& h_in ) const{
uint64_t h = h_in;
return h;
}

int np02_shape_owner::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const{
int err_cnt = 0;
const np02_shape * const shape = get_shape();
if( NULL != shape ){
    if( this != shape->get_shape_owner() ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shape_owner: this=%x != (shape=%x)->owner=%x\n",
             this, shape, shape->get_shape_owner()); 
        }

    const np02_shp_alloc * const shp_alloc = get_shp_alloc();
    const np02_shp_alloc * const shp_shp_alloc = shape->get_shp_alloc();
    if( shp_alloc != shp_shp_alloc ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shape_owner: (this=%x) shp_alloc=%x != (shape=%x)->shp_alloc=%x\n",
             this, shp_alloc, shape, shp_shp_alloc); 
        }
    }
return err_cnt;
}


np02_shape_handle::np02_shape_handle(): np02_shape_owner(), m_shp_alloc(NULL),
    m_alloc_idx(0), m_hndl_type(NP02_SHAPE_HANDLE_TYPE_COUNT),
    m_owner_idx(NP02_SHP_HNDL_OWNER_INVALID_IDX), m_object_ptr(){
}

np02_shape_handle::~np02_shape_handle(){
destruct();
}

void np02_shape_handle::destruct(){
clear_shp_rgn();
}

np02_shp_alloc *np02_shape_handle::get_shp_alloc() const{
return m_shp_alloc;
}

void np02_shape_handle::set_shape(np02_shape *shape){
if((NP02_SHAPE_HANDLE_TYPE_SHAPE != m_hndl_type) ||
    (m_object_ptr.m_shape != shape)){
    clear_shp_rgn();
    m_hndl_type = NP02_SHAPE_HANDLE_TYPE_SHAPE;
    m_object_ptr.m_shape = shape;
    }
AUTO_ASSERT(m_object_ptr.m_shape == shape);
AUTO_ASSERT( ( NULL == shape ) || ( shape->get_shp_alloc() == m_shp_alloc ) );
}

void np02_shape_handle::set_region(np02_region *region){
if((NP02_SHAPE_HANDLE_TYPE_REGION != m_hndl_type) ||
    (m_object_ptr.m_region != region)){
    clear_shp_rgn();
    m_hndl_type = NP02_SHAPE_HANDLE_TYPE_REGION;
    m_object_ptr.m_region = region;
    }
AUTO_ASSERT(m_object_ptr.m_region == region);
AUTO_ASSERT((NULL == region) || (region->get_shp_alloc() == m_shp_alloc));
}

np02_shape *np02_shape_handle::get_shape() const{
np02_shape *shape = NULL;
switch ( m_hndl_type ){
    case NP02_SHAPE_HANDLE_TYPE_SHAPE:
        shape = m_object_ptr.m_shape;
        break;
    case NP02_SHAPE_HANDLE_TYPE_REGION:
        shape = NULL;
        break;
    case NP02_SHAPE_HANDLE_TYPE_COUNT:
    default:
        shape = NULL;
        break;
    }
AUTO_ASSERT( ( NULL == shape ) || ( shape->get_shp_alloc() == m_shp_alloc ) );
return shape;
}

np02_region *np02_shape_handle::get_region() const{
np02_region *region = NULL;
switch ( m_hndl_type ){
    case NP02_SHAPE_HANDLE_TYPE_SHAPE:
        region = NULL;
        break;
    case NP02_SHAPE_HANDLE_TYPE_REGION:
        region = m_object_ptr.m_region;
        break;
    case NP02_SHAPE_HANDLE_TYPE_COUNT:
    default:
        region = NULL;
        break;
    }
AUTO_ASSERT((NULL == region) || (region->get_shp_alloc() == m_shp_alloc));
return region;
}

uint64_t np02_shape_handle::hash( const uint64_t& h_in ) const{
uint64_t h = h_in;
#if defined( CF01_SUPPORT )
h = cf01_obj_hash( h, m_alloc_idx );
h = cf01_obj_hash( h, m_hndl_type );
#else
h += reinterpret_cast<cf01_uint64>(m_alloc_idx);
h ^= ((h << 27) | (h >> 37));
h += reinterpret_cast<cf01_uint64>(m_hndl_type);
h ^= ((h << 27) | (h >> 37));
#endif
switch(m_hndl_type){
    case NP02_SHAPE_HANDLE_TYPE_SHAPE:
        if ( NULL != m_object_ptr.m_shape ){
            h = m_object_ptr.m_shape->hash(h);
            }
        break;
    case NP02_SHAPE_HANDLE_TYPE_REGION:
        if ( NULL != m_object_ptr.m_region ){
            h = m_object_ptr.m_region->hash(h);
            }
        break;
    case NP02_SHAPE_HANDLE_TYPE_FREE_CHAIN: /* fall through */
    case NP02_SHAPE_HANDLE_TYPE_COUNT: /* fall through */
    default:
        break;
    }
return h;
}

int np02_shape_handle::verify_data( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;
err_cnt += np02_shape_owner::verify_data(err_msg,err_msg_capacity,err_msg_pos);
if( NULL == m_shp_alloc ){
    if( 0 != m_alloc_idx ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shape_handle: this=%x  shp_alloc NULL  m_alloc_idx=%i\n",
            this, m_alloc_idx); 
        }
    } 
else{
    const np02_shape_handle *sa_shp_hndl =
        m_shp_alloc->alloc_get_shape_handle_by_idx(m_alloc_idx);
    if( this != sa_shp_hndl ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shape_handle: this=%x  !=  (shp_alloc=%x)->"
            "alloc_get_shape_handle_by_idx(m_alloc_idx=%i)=%x\n",
            this, m_shp_alloc, m_alloc_idx, sa_shp_hndl); 
        }
    }

if( NP02_SHAPE_HANDLE_TYPE_SHAPE == m_hndl_type ){
    if ( NULL != m_object_ptr.m_shape ){
        err_cnt += m_object_ptr.m_shape->verify_data(err_msg,err_msg_capacity,
            err_msg_pos);
        const np02_shape_owner *shp_owner =
            m_object_ptr.m_shape->get_shape_owner();
        if( this != shp_owner ){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shape_handle: this=%x  !=  (shape=%x)->owner=%x\n",
                this, m_object_ptr.m_shape, shp_owner);
            }
        const np02_shp_alloc *shp_shp_alloc =
            m_object_ptr.m_shape->get_shp_alloc();
        if( shp_shp_alloc != m_shp_alloc ){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shape_handle: (this=%x) shp_alloc=%x  !=  "
                "(shape=%x)->shp_alloc=%x\n",
                this, m_shp_alloc, m_object_ptr.m_shape, shp_shp_alloc);
            }
        }
    }
else if( NP02_SHAPE_HANDLE_TYPE_REGION == m_hndl_type ){
    if ( NULL != m_object_ptr.m_region ){
        err_cnt += m_object_ptr.m_region->verify_data(err_msg,err_msg_capacity,
            err_msg_pos);
        const np02_shape_handle *shp_hndl = m_object_ptr.m_region->get_owner();
        if( this != shp_hndl ){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shape_handle: this=%x  !=  (region=%x)->owner=%x\n",
                this, m_object_ptr.m_region, shp_hndl);
            }
        const np02_shp_alloc *rgn_shp_alloc =
            m_object_ptr.m_region->get_shp_alloc();
        if( rgn_shp_alloc != m_shp_alloc ){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shape_handle: (this=%x) shp_alloc=%x  !=  "
                "(region=%x)->shp_alloc=%x\n",
                this, m_shp_alloc, m_object_ptr.m_region, rgn_shp_alloc);
            }
        }
    }
else if( NP02_SHAPE_HANDLE_TYPE_FREE_CHAIN == m_hndl_type ){
    if ( NULL != m_object_ptr.m_free_chain_next ){
        /*  */
        }
    }
else{
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "shape_handle: (this=%x)->m_hndl_type=%i\n", this, m_hndl_type);
    }
return err_cnt;
}

/* free memory */
void np02_shape_handle::clear_shp_rgn(){
switch ( m_hndl_type ){
    case NP02_SHAPE_HANDLE_TYPE_SHAPE:
        if( NULL != m_object_ptr.m_shape ){
            AUTO_ASSERT( m_object_ptr.m_shape->get_shp_alloc() == m_shp_alloc );
            if(NULL == m_shp_alloc){
                delete  m_object_ptr.m_shape;
                }
            else{
                m_shp_alloc->free_shape(m_object_ptr.m_shape);
                }
            m_object_ptr.m_shape = NULL;
            }
        break;
    case NP02_SHAPE_HANDLE_TYPE_REGION:
        if( NULL != m_object_ptr.m_region ){
            AUTO_ASSERT( m_object_ptr.m_region->get_shp_alloc() == m_shp_alloc );
            if(NULL == m_shp_alloc){
                delete  m_object_ptr.m_region;
                }
            else{
                m_shp_alloc->free_region(m_object_ptr.m_region);
                }
            m_object_ptr.m_region = NULL;
            }
        break;
    case NP02_SHAPE_HANDLE_TYPE_COUNT:
    default:
        break;
    }
m_hndl_type = NP02_SHAPE_HANDLE_TYPE_COUNT;
}



int np02_dist_from_xy_xy::verify_result( const np02_shape *a,
    const np02_shape *b, char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const{
AA_INCR_CALL_DEPTH();
int err_cnt = 0;

if(NULL == a ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }
if(NULL == b ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }

if((NULL != a) && (NULL != b)){
    CF01_HASH_CONSISTENCY_CHECK( a->hash(b->hash()) );

    /* bounding box */
    np02_xy bb_a_xy_min, bb_a_xy_max;
    np02_xy bb_b_xy_min, bb_b_xy_max;
    np02_xy bb_ab_xy_min, bb_ab_xy_max;
    a->get_bb(&bb_a_xy_min, &bb_a_xy_max);
    b->get_bb(&bb_b_xy_min, &bb_b_xy_max);
    bb_ab_xy_min.set_x( ( bb_a_xy_min.get_x() < bb_b_xy_min.get_x() ) ?
                          bb_a_xy_min.get_x() : bb_b_xy_min.get_x() );
    bb_ab_xy_min.set_y( ( bb_a_xy_min.get_y() < bb_b_xy_min.get_y() ) ?
                          bb_a_xy_min.get_y() : bb_b_xy_min.get_y() );
    bb_ab_xy_max.set_x( ( bb_a_xy_max.get_x() > bb_b_xy_max.get_x() ) ?
                          bb_a_xy_max.get_x() : bb_b_xy_max.get_x() );
    bb_ab_xy_max.set_y( ( bb_a_xy_max.get_y() > bb_b_xy_max.get_y() ) ?
                          bb_a_xy_max.get_y() : bb_b_xy_max.get_y() );
    
    /* bounding box overlap 
           AAAAAAAAAA            
           A        A            
           A        A            
           A        A            
           A    BBBBBBBBBBBBBBB        
           A    B   A         B  
           A    B   A         B  
           AAAAABAAAA         B  
                BBBBBBBBBBBBBBB 
    */
    const bool bb_overlap = np02_shape::is_bb_overlap( bb_a_xy_min, bb_a_xy_max,
        bb_b_xy_min, bb_b_xy_max );
    
    /* max err distance */
    const double bb_ab_semi_perimeter = 
        ( bb_ab_xy_max.get_x() - bb_ab_xy_min.get_x() ) +
        ( bb_ab_xy_max.get_y() - bb_ab_xy_min.get_y() );
    const double max_err_d = bb_ab_semi_perimeter * np02_shape::m_little_ratio;

    /* distance from near_xy to shape A */
    const double xy_near_a_a_d = a->get_distance_from_xy(m_near_xy);

    if( ( NP02_ANSWER_QUALITY_SMALL_ERROR == m_answer_quality  ) &&
        m_near_xy_defined && ( xy_near_a_a_d > max_err_d ) ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s"
            "  xy_near_a_a_d=%g > max_err_d=%g"
            "  distance_from=%g  near_xy_defined=%i  other_near_xy_defined=%i"
            "  near_xy=(%g,%g)  other_near_xy=(%g,%g)\n", 
            __FILE__, __LINE__, __FUNCTION__,
            xy_near_a_a_d, max_err_d,
            m_distance_from, static_cast<int>(m_near_xy_defined),
            static_cast<int>(m_other_near_xy_defined),
            m_near_xy.get_x(), m_near_xy.get_y(), 
            m_other_near_xy.get_x(), m_other_near_xy.get_y());
        }

    /* distance from near_xy to shape B */
    const double xy_near_b_b_d = b->get_distance_from_xy(m_other_near_xy);

    if( ( NP02_ANSWER_QUALITY_SMALL_ERROR == m_answer_quality  ) &&
        m_other_near_xy_defined && ( xy_near_b_b_d > max_err_d ) ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s"
            "  xy_near_b_b_d=%g > max_err_d=%g"
            "  distance_from=%g  near_xy_defined=%i  other_near_xy_defined=%i"
            "  near_xy=(%g,%g)  other_near_xy=(%g,%g)\n", 
            __FILE__, __LINE__, __FUNCTION__,
            xy_near_b_b_d, max_err_d,
            m_distance_from, static_cast<int>(m_near_xy_defined),
            static_cast<int>(m_other_near_xy_defined),
            m_near_xy.get_x(), m_near_xy.get_y(), 
            m_other_near_xy.get_x(), m_other_near_xy.get_y());
        }

    /* distance from near_xy to shape B */
    const double near_pt_d = m_near_xy.get_distance_to(m_other_near_xy);

    if( bb_overlap ){
        
        }
    else{
        if( ( NP02_ANSWER_QUALITY_SMALL_ERROR == m_answer_quality  ) &&
            ( m_distance_from < -max_err_d ) ){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "data error: %s[%i] %s"
                "  m_distance_from=%g < -max_err_d=%g"
                "  near_xy_defined=%i  other_near_xy_defined=%i"
                "  near_xy=(%g,%g)  other_near_xy=(%g,%g)\n", 
                __FILE__, __LINE__, __FUNCTION__,
                m_distance_from, -max_err_d,
                static_cast<int>(m_near_xy_defined),
                static_cast<int>(m_other_near_xy_defined),
                m_near_xy.get_x(), m_near_xy.get_y(), 
                m_other_near_xy.get_x(), m_other_near_xy.get_y());
            }


        }
    }

AA_DECR_CALL_DEPTH();
return err_cnt;
}


void np02_shape::cos_sin_rot_deg( const double& rot_deg,
    double *cos_rot, double *sin_rot ){
AA_ALWAYS_ASSERT( NULL != cos_rot );
AA_ALWAYS_ASSERT( NULL != sin_rot );
*cos_rot = 1.0;
*sin_rot = 0.0;
if((-360.0 == rot_deg)||(0.0 == rot_deg)||(360.0 == rot_deg)){ 
    *cos_rot = 1.0; *sin_rot = 0.0;
}
else if( (90.0 == rot_deg) || (-270.0 == rot_deg) ){
    *cos_rot = 0.0; *sin_rot = 1.0;
}
else if( (-90.0 == rot_deg) || (270.0 == rot_deg) ){
    *cos_rot = 0.0; *sin_rot = -1.0;
}
else if( (-180.0 == rot_deg) || (180.0 == rot_deg) ){
    *cos_rot = -1.0; *sin_rot = 0.0;
}
else{
    const double rot_rad =
        rot_deg * (3.1415926535897932384626433832795 / 180.0);
    *cos_rot = cos(rot_rad);
    *sin_rot = sin(rot_rad);
    }
AA_XDBG_ASSERT( fabs( *cos_rot ) <= 1.0, AA_DEBUG_LEVEL_2);
AA_XDBG_ASSERT( fabs( *sin_rot ) <= 1.0, AA_DEBUG_LEVEL_2);
AA_XDBG_ASSERT( fabs( 1.0 - ( ( (*cos_rot) * (*cos_rot) ) + 
    ( (*sin_rot) * (*sin_rot) ) ) ) < m_little_ratio_sq, AA_DEBUG_LEVEL_2);
}

np02_shape::np02_shape():m_shape_type(NP02_SHAPE_TYPE_SHAPE),
    m_orientation(NP02_SHAPE_ORIENTATION_COUNT),
    m_invert_status(NP02_SHAPE_INVERT_STATUS_COUNT),
    m_pad(0), m_shape_idx(0), m_shape_owner(NULL),
    m_shp_alloc(NULL),  m_head_loc_grid_node(NULL), 
    m_free_chain_next(NULL){}

double np02_shape::get_small_distance() const{
static const double default_small_distance = np02_shape::m_small_ratio;
const np02_loc_grid *loc_grid = get_loc_grid();
const double small_d = ( NULL == loc_grid ) ?
    default_small_distance : loc_grid->get_small_distance();
return small_d;
}

/* free resources */
void np02_shape::destruct(){
if(NULL != m_head_loc_grid_node ){
    np02_loc_grid *loc_grid = m_head_loc_grid_node->get_loc_grid();
    if(NULL != loc_grid){
        loc_grid->remove_shape_from_loc_grid(this);
        }
    m_head_loc_grid_node = NULL;
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
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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
AA_INCR_CALL_DEPTH();
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
    CF01_HASH_CONSISTENCY_CHECK( shape_a->hash(shape_b->hash()) );

    /* distance between near points */
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
           AAAAAAAAAA            
           A        A            
           A        A            
           A        A            
           A    BBBBBBBBBBBBBBB        
           A    B   A         B  
           A    B   A         B  
           AAAAABAAAA         B  
                BBBBBBBBBBBBBBB 
    */
    const bool bb_overlap = is_bb_overlap( bb_a_xy_min, bb_a_xy_max,
        bb_b_xy_min, bb_b_xy_max );
    
    /* max err distance */
    const double bb_ab_semi_perimeter = 
        ( bb_ab_xy_max.get_x() - bb_ab_xy_min.get_x() ) +
        ( bb_ab_xy_max.get_y() - bb_ab_xy_min.get_y() );
    const double max_err_d = bb_ab_semi_perimeter * 8 * np02_shape::m_little_ratio;

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
            "  distance_from=%f > %f  &&  xy_near_a_b_d_err=%f > %f\n", 
            __FILE__, __LINE__, __FUNCTION__,
            shape_a->get_shape_type(), shape_b->get_shape_type(),
            distance_from, max_err_d, xy_near_a_b_d_err, max_err_d );
        }

    if( (!bb_overlap) && (distance_from < -max_err_d) ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s  a_type=%i  b_type=%i"
            "  (!bb_overlap) && distance_from=%f < %f\n", 
            __FILE__, __LINE__, __FUNCTION__,
            shape_a->get_shape_type(), shape_b->get_shape_type(),
            distance_from, -max_err_d );
        }

    /*   */
    if( d_a_a > max_err_d ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s  a_type=%i  b_type=%i"
            "  d_a_a=%f > %f\n", 
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
            "  distance_from=%f > %f  &&  d_a_b=%f < %f\n", 
            __FILE__, __LINE__, __FUNCTION__,
            shape_a->get_shape_type(), shape_b->get_shape_type(),
            distance_from, max_err_d, d_a_b, -max_err_d );
        }

    if( d_b_b > max_err_d ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s  a_type=%i  b_type=%i"
            "  d_b_b=%f > %f\n", 
            __FILE__, __LINE__, __FUNCTION__,
            shape_a->get_shape_type(), shape_b->get_shape_type(),
            d_b_b, max_err_d );
        }
    /* if( (distance_from < -max_err_d) && (d_b_a > max_err_d) ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s  a_type=%i  b_type=%i"
            "  distance_from=%f < %f  &&  d_b_a=%f > %f\n", 
            __FILE__, __LINE__, __FUNCTION__,
            shape_a->get_shape_type(), shape_b->get_shape_type(),
            distance_from, -max_err_d, d_b_a, max_err_d );
        } */
    if( (distance_from > max_err_d) && (d_b_a < -max_err_d) ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s  a_type=%i  b_type=%i"
            "  distance_from=%f > %f  &&  d_b_a=%f < %f\n", 
            __FILE__, __LINE__, __FUNCTION__,
            shape_a->get_shape_type(), shape_b->get_shape_type(),
            distance_from, max_err_d, d_b_a, -max_err_d );
        }

    if( err_cnt > 0 ){
        std::ostringstream debug_ss;
        debug_ss << "SHAPE_A:\n";
        shape_a->ostream_output(debug_ss);
        debug_ss << "\nSHAPE_B:\n";
        shape_b->ostream_output(debug_ss);
        debug_ss << "\n";
        const std::string debug_ss_str = debug_ss.str();
        np02_snprintf( err_msg, err_msg_capacity, err_msg_pos,
            "%sxy_near_a=[%g,%g]  xy_near_b=[%g,%g]  distance_from=%g\n", 
            debug_ss_str.c_str(), xy_near_a.get_x(), xy_near_a.get_y(),
            xy_near_b.get_x(), xy_near_b.get_y(), distance_from );
        }
    }
AA_DECR_CALL_DEPTH();

return err_cnt;
}

void np02_shape::bb_merge(
    const np02_xy& bb_a_xy_min, const np02_xy& bb_a_xy_max,
    const np02_xy& bb_b_xy_min, const np02_xy& bb_b_xy_max,
    np02_xy *bb_c_xy_min, np02_xy *bb_c_xy_max ){
AA_ALWAYS_ASSERT(NULL != bb_c_xy_min);
AA_ALWAYS_ASSERT(NULL != bb_c_xy_max);
bb_c_xy_min->set_x( (bb_a_xy_min.get_x() < bb_b_xy_min.get_x()) ? 
                     bb_a_xy_min.get_x() : bb_b_xy_min.get_x() );
bb_c_xy_min->set_y( (bb_a_xy_min.get_y() < bb_b_xy_min.get_y()) ? 
                     bb_a_xy_min.get_y() : bb_b_xy_min.get_y() );
bb_c_xy_max->set_x( (bb_a_xy_max.get_x() > bb_b_xy_max.get_x()) ? 
                     bb_a_xy_max.get_x() : bb_b_xy_max.get_x() );
bb_c_xy_max->set_y( (bb_a_xy_max.get_y() > bb_b_xy_max.get_y()) ? 
                     bb_a_xy_max.get_y() : bb_b_xy_max.get_y() );
}

/*
>0  ==> bounding boxes separated
<0  ==> bounding boxes overlap

*/
double np02_shape::get_bb_sep_dist_manhattan(
    const np02_xy& bb_a_xy_min, const np02_xy& bb_a_xy_max,
    const np02_xy& bb_b_xy_min, const np02_xy& bb_b_xy_max ){
double s = 0.0;

AA_ALWAYS_ASSERT(false); /* not implemented */

return s;
}

uint64_t np02_shape::hash( const uint64_t& h_in ) const{
uint64_t h = h_in;
#if defined( CF01_SUPPORT )
h = cf01_obj_hash( h, m_shape_type );
h = cf01_obj_hash( h, m_orientation );
h = cf01_obj_hash( h, m_invert_status );
AUTO_ASSERT( 0 == m_pad );
h = cf01_obj_hash( h, m_shape_idx );
if( NULL != m_shape_owner ){
    h = cf01_obj_hash( h, static_cast<uint8_t>(1) );
    }
if( NULL != m_shp_alloc ){
    h = cf01_obj_hash( h, static_cast<uint8_t>(1) );
    }
#else
h += reinterpret_cast<cf01_uint64>(m_shape_type);
h ^= ((h << 27) | (h >> 37));
h += reinterpret_cast<cf01_uint64>(m_orientation);
h ^= ((h << 27) | (h >> 37));
h += reinterpret_cast<cf01_uint64>(m_invert_status);
h ^= ((h << 27) | (h >> 37));
AUTO_ASSERT( 0 == m_pad );
h += reinterpret_cast<cf01_uint64>(m_shape_idx);
h ^= ((h << 37) | (h >> 27));
if( NULL != m_shape_owner ){
    h += 1;
    h ^= ((h << 47) | (h >> 17));
    }
if( NULL != m_shp_alloc ){
    h += 1;
    h ^= ((h << 47) | (h >> 17));
    }
#endif

if( NULL != m_head_loc_grid_node ){ 
    h = m_head_loc_grid_node->hash( h );
    const np02_loc_grid_node *p = m_head_loc_grid_node->get_s_prev();
    while( NULL != p ){ 
        h = p->hash( h ); 
        p = p->get_s_prev();
        }
    const np02_loc_grid_node *n = m_head_loc_grid_node->get_s_next();
    while( NULL != p ){ 
        h = n->hash( h ); 
        n = n->get_s_next();
        }
    }
return h;
}

int np02_shape::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const{
int err_cnt = 0;
if( m_shape_type > NP02_SHAPE_TYPE_COUNT){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shape: this=%x  m_shape_type=%i\n", this, m_shape_type); 
    }

if( m_orientation > NP02_SHAPE_ORIENTATION_COUNT ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shape: this=%x  m_orientation=%i\n", this, m_orientation); 
    }

if( m_invert_status > NP02_SHAPE_INVERT_STATUS_COUNT ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shape: this=%x  m_invert_status=%i\n", this, m_invert_status); 
    }

if( 0 != m_pad ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shape: this=%x  m_pad=%i\n", this, m_pad); 
    }

if(m_shape_idx > NP02_SHAPE_MAX_IDX){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shape: this=%x  m_shape_idx=%i\n", this, m_shape_idx); 
    }

if( NULL != m_shape_owner ){
    const np02_shape * const ownr_shp = m_shape_owner->get_shape();
    if( this != ownr_shp ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shape: this=%x != (owner=%x)->shape=%x\n", this,
            m_shape_owner, ownr_shp); 
        }
    const np02_shp_alloc * const ownr_shp_alloc=m_shape_owner->get_shp_alloc();
    if( m_shp_alloc != ownr_shp_alloc ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shape: this=%x shp_alloc=%x  !=  (owner=%x)->shp_alloc=%x\n",
            this, m_shp_alloc, m_shape_owner, ownr_shp_alloc); 
        }
    }

if( NULL != m_head_loc_grid_node){
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
os << "<shape_type>" << m_shape_type << "</shape_type>\n";
os << "<orientation>" << m_orientation << "</orientation>\n";
os << "<invert_status>" << m_invert_status << "</invert_status>\n";
AUTO_ASSERT( 0 == m_pad );
os << "<shape_idx>" << m_shape_idx << "</shape_idx>\n";
os << std::hex;
os << "<shape_owner>" << m_shape_owner << "</shape_owner>\n";
os << "<shp_alloc>" << m_shp_alloc << "</shp_alloc>\n";
os << "<head_loc_grid_node>" << m_head_loc_grid_node << "</head_loc_grid_node>\n";
os << "<m_free_chain_next>" << m_free_chain_next << "</m_free_chain_next>\n";
os << std::dec;
os << "</shape>\n";
return os;
}

void np02_shape::write_bmp_file(const np02_xy& xy_min,
    const double& pixel_num, const np02_bmp_color& color,
    np02_bmp_file *bmp_file) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AA_ALWAYS_ASSERT(NULL != bmp_file);
const int32_t& width_px = bmp_file->get_width_px();
const int32_t& height_px = bmp_file->get_height_px();

np02_xy bb_xy_min, bb_xy_max;
get_bb(&bb_xy_min, &bb_xy_max);

const int32_t i_min_a = static_cast<int32_t>(
    (bb_xy_min.get_x() - xy_min.get_x()) * pixel_num);
const int32_t i_min = (i_min_a < 0 ) ? 0 : 
    ( i_min_a >= width_px ) ? (width_px-1) : i_min_a;

const int32_t i_max_a = static_cast<int32_t>(
    (bb_xy_max.get_x() - xy_min.get_x()) * pixel_num);
const int32_t i_max = (i_max_a < 0 ) ? 0 : 
    ( i_max_a >= width_px ) ? (width_px-1) : i_max_a;

const int32_t j_min_a = static_cast<int32_t>(
    (bb_xy_min.get_y() - xy_min.get_y()) * pixel_num);
const int32_t j_min = (j_min_a < 0 ) ? 0 : 
    ( j_min_a >= height_px ) ? (height_px-1) : j_min_a;

const int32_t j_max_a = static_cast<int32_t>(
    (bb_xy_max.get_y() - xy_min.get_y()) * pixel_num);
const int32_t j_max = (j_max_a < 0 ) ? 0 : 
    ( j_max_a >= height_px ) ? (height_px-1) : j_max_a;

np02_xy xy;
for(int32_t i = i_min; i <= i_max; ++i){
    AA_INCR_CALL_DEPTH();
    xy.set_x(xy_min.get_x() + (static_cast<double>(i) / pixel_num));
    for(int32_t j = j_min; j <= j_max; ++j){
        AA_INCR_CALL_DEPTH();
        CF01_HASH_CONSISTENCY_CHECK( this->hash() );
        xy.set_y(xy_min.get_y() + (static_cast<double>(j) / pixel_num));
        const double d = get_distance_from_xy(xy);
        if(d <= 0.0 ){
            bmp_file->draw_pixel(i,j,color);
            }
        AA_DECR_CALL_DEPTH();
        }
    AA_DECR_CALL_DEPTH();
    }
AA_DECR_CALL_DEPTH();
}

void np02_shape::write_dxf_file(const std::string& layer,
    const uint8_t& color, np02_dxf_file *dxf_file) const{

}


uint64_t np02_shape_vec_hash( const np02_shape_vec& v,
    const uint64_t& h_in ){
uint64_t h = h_in;
for( np02_shape_vec_citr v_itr = v.begin(); v_itr != v.end(); ++v_itr ){
    const np02_shape *s = *v_itr;
    if ( NULL == s) {
        static const uint8_t zero = 0;
        h = cf01_obj_hash( h, zero );
        }
    else{
        h = s->hash(h);
        }
    }
return h;
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
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
const double dx = xy.get_x() - m_ctr.get_x();
const double dy = xy.get_y() - m_ctr.get_y();
const double dsq = (dx*dx) + (dy*dy);
const double d_to_ctr = sqrt(dsq);
const double d = d_to_ctr - m_radius;
if(NULL != near_xy){
    if(m_radius < np02_shape::m_teensy_ratio) {
        *near_xy = m_ctr; }
    else if(d_to_ctr < np02_shape::m_teensy_ratio) {
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
AA_DECR_CALL_DEPTH();
return d;
}

double np02_circle::get_distance_from_line_seg_ab(const np02_xy& xy_a,
    const np02_xy& xy_b, np02_xy *near_xy,
    np02_xy *other_near_xy) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AUTO_ASSERT(m_radius >= 0.0);
double dx, dy, dsq, d;
double d_ctr = 0.0;

/* forward unit vector A->B */
const double dx_ab = xy_b.get_x() - xy_a.get_x();
const double dy_ab = xy_b.get_y() - xy_a.get_y();
const double dsq_ab = (dx_ab * dx_ab) + (dy_ab * dy_ab);
np02_xy fwd_ab(1.0,0.0);
if(dsq_ab > np02_shape::m_teensy_ratio){
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
        if(m_radius < np02_shape::m_teensy_ratio) {
            *near_xy = m_ctr; }
        else if(d_ctr < np02_shape::m_teensy_ratio) {
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
        if(m_radius < np02_shape::m_teensy_ratio) {
            *near_xy = m_ctr; }
        else if(d_ctr < np02_shape::m_teensy_ratio) {
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
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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
    if( d_ctr < np02_shape::m_teensy_ratio){
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

double np02_circle::get_distance_from_arc(const np02_arc *a,
    np02_xy *near_xy, np02_xy *arc_near_xy) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AA_ALWAYS_ASSERT(NULL != a);
double d = a->get_distance_from_circle( this, arc_near_xy, near_xy);
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
        AUTO_ASSERT(fabs(dd - fabs(d_ctr)) <= 
            fabs((dd + fabs(d_ctr))*np02_shape::m_little_ratio) );
        np02_xy unit_v(1.0,0.0);
        if( dd < np02_shape::m_teensy_ratio){
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
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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
        AUTO_ASSERT(fabs(dd - fabs(d_ctr)) <=
            fabs((dd + fabs(d_ctr))*np02_shape::m_little_ratio) );
        np02_xy unit_v(1.0,0.0);
        if( dd < np02_shape::m_teensy_ratio){
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
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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
        AUTO_ASSERT(fabs(dd - d_ctr) <=
            fabs((dd + d_ctr)*np02_shape::m_little_ratio) );
        np02_xy unit_v(1.0,0.0);
        if( dd < np02_shape::m_teensy_ratio){
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
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
if(NULL != s){
    const np02_circle *c = dynamic_cast<const np02_circle *>(s);
    const np02_arc *a = dynamic_cast<const np02_arc *>(s);
    const np02_line_seg *n = dynamic_cast<const np02_line_seg *>(s);
    const np02_rect *r = dynamic_cast<const np02_rect *>(s);
    const np02_polygon *p = dynamic_cast<const np02_polygon *>(s);
    const np02_spline *i = dynamic_cast<const np02_spline *>(s);
    if(NULL != c){ d = get_distance_from_circle(c, near_xy, other_near_xy); }
    else if(NULL!=a){ d=get_distance_from_arc(a, near_xy, other_near_xy);}
    else if(NULL!=n){ d=get_distance_from_line_seg(n, near_xy, other_near_xy);}
    else if(NULL!=r){ d=get_distance_from_rect(r, near_xy, other_near_xy); }
    else if(NULL!=p){ d=get_distance_from_polygon(p, near_xy, other_near_xy); }
    else if(NULL!=i){ AA_ALWAYS_ASSERT(false); /* spline query not supported */ }
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
if( ( 0.0 != rot_deg ) && ( -360.0 != rot_deg ) && ( 360.0 != rot_deg ) ){
    const np02_xy rot_arm_initial(m_ctr.get_x() - rot_ctr.get_x(),
        m_ctr.get_y() - rot_ctr.get_y());
    if( (rot_arm_initial.get_x() != 0.0) || (rot_arm_initial.get_y() != 0.0) ){
        double cos_rot = 1.0;
        double sin_rot = 0.0;
        cos_sin_rot_deg( rot_deg, &cos_rot, &sin_rot );
        const np02_xy rot_arm_final(
            (rot_arm_initial.get_x() * cos_rot) -
            (rot_arm_initial.get_y() * sin_rot),
            (rot_arm_initial.get_x() * sin_rot) +
            (rot_arm_initial.get_y() * cos_rot));
        m_ctr.set_x(rot_ctr.get_x() + rot_arm_final.get_x());
        m_ctr.set_y(rot_ctr.get_y() + rot_arm_final.get_y());
        }
    }
}

uint64_t np02_circle::hash( const uint64_t& h_in ) const{
uint64_t h = np02_shape::hash( h_in );
h = m_ctr.hash( h );
#if defined( CF01_SUPPORT )
h = cf01_obj_hash( h, m_radius );
#else
h += reinterpret_cast<cf01_uint64>(m_radius);
h ^= ((h << 27) | (h >> 37));
#endif
return h;
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
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AUTO_ASSERT(NULL != dxf_file);
dxf_file->draw_circle( layer, m_ctr.get_x(), m_ctr.get_y(), m_radius, color );
AA_DECR_CALL_DEPTH();
}

size_t np02_circle::circle_circle_intersect(const np02_circle *c,
    np02_xy *xy_intsct_0, np02_xy *xy_intsct_1) const{
size_t intsct_count = 0;
if( NULL != c ){
    intsct_count = circle_circle_intsct( m_ctr, m_radius, c->m_ctr,
        c->m_radius, xy_intsct_0, xy_intsct_1);
    }
return intsct_count;
}

/* return number of intersections */
size_t np02_circle::circle_circle_intsct(const np02_xy& ctr_a,
    const double& r_a, const np02_xy& ctr_b, const double& r_b,
    np02_xy *xy_intsct_0, np02_xy *xy_intsct_1){
AA_INCR_CALL_DEPTH();
AUTO_ASSERT(r_a >= 0.0);
AUTO_ASSERT(r_b >= 0.0);
size_t intsct_count = 0;
if( ctr_a == ctr_b ){
    if( r_a == r_b ){
        if( 0.0 == r_a ){
            intsct_count = 1;
            if( NULL != xy_intsct_0 ){
                *xy_intsct_0 = ctr_a;
                }
            }
        else{
            intsct_count = std::numeric_limits<size_t>::max();
            if( NULL != xy_intsct_0 ){
                xy_intsct_0->set_x( ctr_a.get_x() - r_a );
                xy_intsct_0->set_y( ctr_a.get_y() );
                }
            if( NULL != xy_intsct_1 ){
                xy_intsct_1->set_x( ctr_a.get_x() + r_a );
                xy_intsct_1->set_y( ctr_a.get_y() );
                }
            }
        }
    else{
        intsct_count = 0;
        }    
    }
else{
    const double dx_a_b = ctr_b.get_x() - ctr_a.get_x();
    const double dy_a_b = ctr_b.get_y() - ctr_a.get_y();
    np02_xy lambda_a_b(1.0, 0.0);
    const double d_a_b_sq = (dx_a_b*dx_a_b) + (dy_a_b*dy_a_b);
    double d_a_b = 0.0; 
    if( d_a_b_sq > 0.0 ){
        d_a_b = sqrt(d_a_b_sq);
        lambda_a_b.set_x(dx_a_b/d_a_b);
        lambda_a_b.set_y(dy_a_b/d_a_b);
        }

    if( ( d_a_b > (r_a + r_b) ) ||
        ( r_a > (r_b + d_a_b) ) ||
        ( r_b > (d_a_b + r_a) ) ){
        intsct_count = 0;
        }
    else{
        /*
                g > 0                  g < 0                     g < -d_a_b
                          *                      *          *
                       / /|                   /  |\         |\ \
                    /   / |                /     | \        | \   \
              r_a/     /  |h         r_a/        |h \      h|  \     \r_b
              /    r_b/   |          /           |   \r_b   |   \r_a    \
           /         /    |       /              |  g \     |    \         \
        *-----------*-----*    *-----------------*----*     *-----*-----------*
        A   d_a_b   B  g  D    A                 D    B     D     A   d_a_b   B 
                               + - - - -d_a_b- - - - -+     + - - - -g- - - - +      
    
        h^2 = b^2 - g^2 = a^2 - (d+g)^2
              b^2 - g^2 = a^2 - d^2 - 2dg - g^2
                    2dg = a^2 - b^2 - d^2
                      g = (a^2 - b^2 - d^2)/2d
        */
        const double r_a_sq = r_a * r_a;
        const double r_b_sq = r_b * r_b;
        const double g = (d_a_b>0.0) ?
            ((r_a_sq-(r_b_sq+d_a_b_sq))/(2.0*d_a_b)) : r_b;
        const double g_sq = g * g;
        const double h_sq = r_b_sq - g_sq;
        const double h = ( h_sq > 0.0 ) ? sqrt( h_sq ) : 0.0;
        intsct_count = ( h > 0.0 ) ? 2 : 1;
        const np02_xy d( ctr_b.get_x() + (g * lambda_a_b.get_x()), 
                         ctr_b.get_y() + (g * lambda_a_b.get_y()) );
        const double h_x = h * lambda_a_b.get_x();
        const double h_y = h * lambda_a_b.get_y();
        if( NULL != xy_intsct_0 ){
            xy_intsct_0->set_x(d.get_x() - h_y);
            xy_intsct_0->set_y(d.get_y() + h_x);
            }
        if( ( h > 0 ) && ( NULL != xy_intsct_1 ) ){
            xy_intsct_1->set_x(d.get_x() + h_y);
            xy_intsct_1->set_y(d.get_y() - h_x);
            }
        }
    }

AUTO_ASSERT( 0 == verify_data_circle_circle_intsct_result( ctr_a,
    r_a, ctr_b, r_b, xy_intsct_0, xy_intsct_1, intsct_count,
    AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ) );

AA_DECR_CALL_DEPTH();
return intsct_count;
}

/* return 0 => result is correct */
int np02_circle::verify_data_circle_circle_intsct_result( const np02_xy& ctr_a,
    const double& r_a, const np02_xy& ctr_b, const double& r_b,
    np02_xy *xy_intsct_0, np02_xy *xy_intsct_1, const size_t& intsct_count,
    char *err_msg, const size_t err_msg_capacity, size_t *err_msg_pos ){
int err_cnt = 0;

const double r_a_sq = r_a * r_a;
const double r_b_sq = r_b * r_b;
const double r_ab_sq_sum = r_a_sq + r_b_sq;
const double max_d_err_sq = ( r_ab_sq_sum > 1.0 ) ? 
    ( r_ab_sq_sum * m_little_ratio_sq ) : m_little_ratio_sq;

if( ( intsct_count > 0 ) && ( NULL != xy_intsct_0 ) ){
    const double dsq_a_0 = xy_intsct_0->get_dsq_to( ctr_a );
    const double dsq_a_0_err = fabs( dsq_a_0 - r_a_sq );
    if( dsq_a_0_err > max_d_err_sq ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_a_0=%g != r_a=%g\n", sqrt(dsq_a_0), r_a ); 
        }

    const double dsq_b_0 = xy_intsct_0->get_dsq_to( ctr_b );
    const double dsq_b_0_err = fabs( dsq_b_0 - r_b_sq );
    if( dsq_b_0_err > max_d_err_sq ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_b_0=%g != r_b=%g\n", sqrt(dsq_b_0), r_b ); 
        }
    }

if( ( intsct_count > 1 ) && ( NULL != xy_intsct_1 ) ){
    const double dsq_a_1 = xy_intsct_1->get_dsq_to( ctr_a );
    const double dsq_a_1_err = fabs( dsq_a_1 - r_a_sq );
    if( dsq_a_1_err > max_d_err_sq ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_a_1=%g != r_a=%g\n", sqrt(dsq_a_1), r_a ); 
        }

    const double dsq_b_1 = xy_intsct_1->get_dsq_to( ctr_b );
    const double dsq_b_1_err = fabs( dsq_b_1 - r_b_sq );
    if( dsq_b_1_err > max_d_err_sq ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_b_1=%g != r_b=%g\n", sqrt(dsq_b_1), r_b ); 
        }
    }

if( err_cnt > 0 ){    
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "verify_data_circle_circle_intsct_result( ctr_a=[%g,%g],"
        " r_a=%g, ctr_b=[%g,%g], r_b=%g, xy_intsct_0=[%g,%g],"
        " xy_intsct_1=[%g,%g], intsct_count=%i\n",
        ctr_a.get_x(), ctr_a.get_y(), r_a, ctr_b.get_x(), ctr_b.get_y(), r_b,
        ( NULL == xy_intsct_0 ) ? 0.0 : xy_intsct_0->get_x(),
        ( NULL == xy_intsct_0 ) ? 0.0 : xy_intsct_0->get_y(),
        ( NULL == xy_intsct_1 ) ? 0.0 : xy_intsct_1->get_x(),
        ( NULL == xy_intsct_1 ) ? 0.0 : xy_intsct_1->get_y(),
        intsct_count );
    }

return err_cnt;
}

np02_arc::np02_arc():m_ctr(0.0, 0.0), m_radius(0.0), m_start_angle_deg(0.0),
    m_end_angle_deg(0.0), m_width(0.0), m_p_0(0.0, 0.0), m_p_1(0.0, 0.0),
    m_fwd_0(0.0, 1.0), m_fwd_1(0.0, 1.0), m_fwd_dot_0(0.0), m_fwd_dot_1(0.0),
    m_bb_xy_min(0.0, 0.0), m_bb_xy_max(0.0, 0.0)
{
set_shape_type(NP02_SHAPE_TYPE_ARC);
}

np02_arc::~np02_arc(){}

void np02_arc::init(const init_params& prm){
AA_INCR_CALL_DEPTH();
const np02_xy& ctr = prm.m_ctr; 
const double& radius = prm.m_radius; 
const double& start_angle_deg = prm.m_start_angle_deg; 
const double& end_angle_deg = prm.m_end_angle_deg;
const double& width = prm.m_width;
AUTO_ASSERT(radius >= 0.0);
AUTO_ASSERT(start_angle_deg > -720.0);
AUTO_ASSERT(start_angle_deg < 720.0);
AUTO_ASSERT(end_angle_deg > -720.0);
AUTO_ASSERT(end_angle_deg < 720.0);
AUTO_ASSERT(width >= 0.0);
m_ctr = ctr;
m_radius = radius;
m_start_angle_deg = start_angle_deg;
if( ( m_start_angle_deg <= -360.0 ) || ( m_start_angle_deg >= 360.0 ) ){
    m_start_angle_deg = fmod( m_start_angle_deg, 360.0 );
    }
double angle_sz = end_angle_deg - m_start_angle_deg;
if( ( angle_sz >= 0.0 ) && ( angle_sz < 360.0 ) ){
    m_end_angle_deg = end_angle_deg;
    }
else{
    angle_sz = fmod( angle_sz, 360.0 );
    if( angle_sz < 0.0 ){
        angle_sz += 360.0;
        }
    m_end_angle_deg = m_start_angle_deg + angle_sz;
    }
AUTO_ASSERT(m_start_angle_deg > -360.0);
AUTO_ASSERT(m_start_angle_deg < 720.0);
AUTO_ASSERT(m_start_angle_deg <= m_end_angle_deg);
AUTO_ASSERT(m_end_angle_deg > -360.0);
AUTO_ASSERT(m_end_angle_deg < 720.0);
AUTO_ASSERT(m_end_angle_deg < (m_start_angle_deg + 360.0));
m_width = width;
double cos_start_angle_deg = 1.0;
double sin_start_angle_deg = 0.0;
cos_sin_rot_deg( m_start_angle_deg, &cos_start_angle_deg,
    &sin_start_angle_deg );
double cos_end_angle_deg = 1.0;
double sin_end_angle_deg = 0.0;
cos_sin_rot_deg( m_end_angle_deg, &cos_end_angle_deg,
    &sin_end_angle_deg );
const double dx0 = radius * cos_start_angle_deg;
const double dy0 = radius * sin_start_angle_deg;
const double dx1 = radius * cos_end_angle_deg;
const double dy1 = radius * sin_end_angle_deg;
m_p_0.set_x( ctr.get_x() + dx0 );
m_p_0.set_y( ctr.get_y() + dy0 );
m_p_1.set_x( ctr.get_x() + dx1 );
m_p_1.set_y( ctr.get_y() + dy1 );
m_fwd_0.set_x( -sin_start_angle_deg );
m_fwd_0.set_y( cos_start_angle_deg );
m_fwd_1.set_x( -sin_end_angle_deg );
m_fwd_1.set_y( cos_end_angle_deg );
m_fwd_dot_0 = m_fwd_0.dot(m_p_0);
m_fwd_dot_1 = m_fwd_1.dot(m_p_1);

/* bounding box */
init_bb();

AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_1);
AA_DECR_CALL_DEPTH();
}

void np02_arc::init_force_p01(const init_params& prm, 
    const init3pt_aux_params& aux_prm) {
AA_INCR_CALL_DEPTH();
init(prm);
m_p_0 = aux_prm.m_p_0;
m_p_1 = aux_prm.m_p_1;
init_bb();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_1);
AA_DECR_CALL_DEPTH();
}

/*
Generate init_params and init3pt_aux_params from 3 points defining an arc.

https://mathworld.wolfram.com/Circumcircle.html

                 C
                 *
              /   \
           /       \
        /           \
    A*---------------*B
*/
void np02_arc::init3pt_to_init_params( const init3pt_params& init3pt_prm,
    init_params *init_prm, init3pt_aux_params *aux_params ){
AA_INCR_CALL_DEPTH();
const np02_xy& pa = init3pt_prm.m_pa; 
const np02_xy& pb = init3pt_prm.m_pb;
const np02_xy& pc = init3pt_prm.m_pc;
const double& width = init3pt_prm.m_width;
const double& max_radius = init3pt_prm.m_max_radius;
AA_ALWAYS_ASSERT(width >= 0.0);
AA_ALWAYS_ASSERT(max_radius > 0.0);

/* x1 = 0 */
/* y1 = 0 */
const double x2 = pb.get_x() - pa.get_x();
const double y2 = pb.get_y() - pa.get_y();
const double x3 = pc.get_x() - pa.get_x();
const double y3 = pc.get_y() - pa.get_y();
const double z3 = (x3*x3) + (y3*y3);

/*
    |  x1   y1   1  |
a = |  x2   y2   1  |
    |  x3   y3   1  |
*/
const double a = (x2*y3) - (x3*y2);

/*
    a = (x2*y3) - (x3*y2) determins CCW vs CW orientation

                  C                          B     
                  * p_1                      *     
               /   \                      /   \    
            /ctr    \                  / ctr   \   
         /    +      \              /     +     \                    2
     A*---------------*B        A*---------------*C        A*--------B------*C          
  p_0       a > 0            p_1        a < 0      p_0             a = 0
        ccwabc = true              ccwabc = false
*/
const bool ccwabc = ( a < 0 ) ? false : true;
const np02_xy& p_0 = ccwabc ? pa : pc;
const np02_xy& p_1 = ccwabc ? pc : pa;

np02_xy ctr(0.0,0.0);
double r = np02_shape::m_googol;
bool is_straight = false;
if (fabs(a) < np02_shape::m_teensy_ratio) {
    is_straight = true;
    }
else{
    /*
           |  x1^2+y1^2   y1   1  |
    b_x = -|  x2^2+y2^2   y2   1  |
           |  x2^2+y2^2   y3   1  |
    */
    const double z2 = (x2*x2) + (y2*y2);
    const double b_x = (z3*y2) - (z2*y3);
    
    /*
           |  x1^2+y1^2   x1   1  |
    b_y =  |  x2^2+y2^2   x2   1  |
           |  x2^2+y2^2   x3   1  |
    */
    const double b_y = (z2*x3) - (z3*x2);
    
    /*
           |  x1^2+y1^2   x1   y1  |
    c  =  -|  x2^2+y2^2   x2   y2  | = 0
           |  x2^2+y2^2   x3   y3  |
    */

    /* circumcenter */
    ctr.set_x( pa.get_x() - b_x/(2.0 * a) );
    ctr.set_y( pa.get_y() - b_y/(2.0 * a) );

    /* radius = sqrt(b_x^2 + b_y^2 -4*a*c)/(2*abs(a))*/
    const double b_x_2 = b_x * b_x;
    const double b_y_2 = b_y * b_y;
    r = sqrt( b_x_2 + b_y_2 )/(2.0 * fabs(a));
    if( r > max_radius ){
        is_straight = true;
        }
    else{
        is_straight = false;
        }
    }

if( is_straight ){
    /*
                                                   * C
                                                   |
    ctr                 ccwabc=true                |
     *---------------------------------------------* mid
                                                   |
                                                   |
                                                   * A
        * C
        |
        |              ccwabc=false
    mid *---------------------------------------------* ctr 
        |
        |
        * A
    */
    const np02_xy mid( (pa.get_x() + pc.get_x() )/2.0,
                       (pa.get_y() + pc.get_y() )/2.0 );
    const double d_ac = sqrt(z3);
    np02_xy fwd_ac(1.0,0.0);
    if (d_ac > 0.0) {
        fwd_ac.set_xy(x3/d_ac, y3/d_ac);
        }
    const np02_xy unit_v_mid_to_ctr(
        ( ccwabc ) ? -fwd_ac.get_y() : fwd_ac.get_y(),
        ( ccwabc ) ? fwd_ac.get_x() : -fwd_ac.get_x() );
    r = max_radius;
    const double r_sq = r * r;
    const double e_sq = r_sq - (z3/4.0);
    const double e = (e_sq > 0.0 ) ? sqrt(e_sq) : 0.0;
    ctr.set_x( mid.get_x() + ( e * unit_v_mid_to_ctr.get_x() ) );
    ctr.set_y( mid.get_y() + ( e * unit_v_mid_to_ctr.get_y() ) );
    }

const np02_xy dxy_ctr_to_p_0( p_0.get_x() - ctr.get_x(), 
                              p_0.get_y() - ctr.get_y());
const np02_xy dxy_ctr_to_p_1( p_1.get_x() - ctr.get_x(), 
                              p_1.get_y() - ctr.get_y() );    
const double start_angle_rad = atan2( dxy_ctr_to_p_0.get_y(), dxy_ctr_to_p_0.get_x() ); 
const double end_angle_rad   = atan2( dxy_ctr_to_p_1.get_y(), dxy_ctr_to_p_1.get_x() );
double start_angle_deg = m_degrees_per_radian * start_angle_rad;
double end_angle_deg = m_degrees_per_radian * end_angle_rad;
if( end_angle_deg < start_angle_deg ){
    if( (start_angle_deg + end_angle_deg) < 0.0 ){
        end_angle_deg += 360.0;
        }
    else{
        start_angle_deg -= 360.0;
        }
    }

if( NULL != init_prm ){
    init_prm->m_ctr = ctr;
    init_prm->m_radius = r;
    init_prm->m_start_angle_deg = start_angle_deg;
    init_prm->m_end_angle_deg = end_angle_deg;
    init_prm->m_width = width;
    }

if( NULL != aux_params ){
    aux_params->m_p_0 = p_0;
    aux_params->m_p_1 = p_1;
    aux_params->m_is_straight = is_straight;
    }

AA_DECR_CALL_DEPTH();
}

bool np02_arc::arc_is_approx_straight_line_seg() const{
static const double almost_straight_angle_deg = 1e-6;
const double angle_range = m_end_angle_deg-m_start_angle_deg; 
return ( angle_range < almost_straight_angle_deg ) ? true : false;
}


/*
                           
                                  true     .
         true                            .
                         ***           .
  .                 *           *    .
      .          *                 +P0
          .    *       true       .      
              +P1               .
                   .           .        false
    false               .    .
                           +ctr
         
                         false
*/
bool np02_arc::is_in_perp_zone(const np02_xy& p) const{
const bool in_perp_0 = !is_in_zone_0(p);
const bool in_perp_1 = !is_in_zone_1(p);
const bool in_perp_zone = ( is_less_than_half_circle() ) ?
    ( in_perp_0 && in_perp_1 ) : ( in_perp_0 || in_perp_1 );
return in_perp_zone;
}

/* result > 0 => p is inside perp zone */
double np02_arc::distance_in_perp_zone(const np02_xy& p) const{
const double d_in_perp_0 = -distance_in_zone_0(p);
const double d_in_perp_1 = -distance_in_zone_1(p);
double d = 0.0;
if( is_less_than_half_circle() ){
    d = ( d_in_perp_0 < d_in_perp_1 ) ? d_in_perp_0 : d_in_perp_1;
    }
else{
    d = ( d_in_perp_0 > d_in_perp_1 ) ? d_in_perp_0 : d_in_perp_1;
    }
return d;
}

void np02_arc::get_bb(np02_xy *xy_min, np02_xy *xy_max) const{
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != xy_min);
AA_ALWAYS_ASSERT(NULL != xy_max);
*xy_min = m_bb_xy_min;
*xy_max = m_bb_xy_max;
AA_DECR_CALL_DEPTH();
}

void np02_arc::get_loc_grid_indices_for_init(
    const np02_loc_grid_dim& loc_grid_dim, const double& extra_search_d,
    np02_uint16_pair_vec *index_vec) const{
AA_INCR_CALL_DEPTH();
AUTO_ASSERT(m_radius >= 0.0);
AUTO_ASSERT(m_width >= 0.0);
AUTO_ASSERT(loc_grid_dim.get_sq_size() > 0.0);
AUTO_ASSERT(extra_search_d >= 0.0);
AA_ALWAYS_ASSERT(NULL != index_vec);

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

np02_xy xy(0.0,0.0);
bool i_on_loc_grid_edge, j_on_loc_grid_edge, should_add;
for( i = ij_min.first; i <= ij_max.first; ++i ){
    xy.set_x(loc_grid_dim.get_sq_ctr_x(i));
    i_on_loc_grid_edge = ((0 == i) || (loc_grid_dim.get_w() == (i+1))) ? 
        true : false;
    for( j = ij_min.second; j <= ij_max.second; ++j ){
        j_on_loc_grid_edge = ((0 == j) || (loc_grid_dim.get_h() == (j+1))) ?
            true : false;
        if(i_on_loc_grid_edge || j_on_loc_grid_edge){
            should_add = true;
            }
        else{
            xy.set_y(loc_grid_dim.get_sq_ctr_y(j));
            const double d_from = get_distance_from_xy(xy, NULL);
            should_add = (d_from <= d_threshold) ? true : false;
            }
        if(should_add){
            index_vec->push_back(np02_uint16_pair(i,j));
            }
        }
    }

AA_DECR_CALL_DEPTH();
}

double np02_arc::get_distance_from_xy(const np02_xy& xy,
    np02_xy *near_xy) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
double d = 0.0;
if( arc_is_approx_straight_line_seg() ){
    np02_line_seg straight_seg;
    straight_seg.init( m_p_0, m_p_1, m_width );
    d = straight_seg.get_distance_from_xy(xy, near_xy );
    }
else{
    const double hw = m_width / 2.0;
    d = get_distance_from_xy_hw(xy, hw, near_xy);
    }
AA_DECR_CALL_DEPTH();
return d;
}

double np02_arc::get_distance_from_line_seg_ab(const np02_xy& xy_a,
    const np02_xy& xy_b, np02_xy *near_xy,
    np02_xy *other_near_xy) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AUTO_ASSERT(m_radius >= 0.0);
double d = 0.0;
np02_line_seg line_seg_ab;
line_seg_ab.init( xy_a, xy_b, 0.0);
if( arc_is_approx_straight_line_seg() ){
    np02_line_seg straight_seg;
    straight_seg.init( m_p_0, m_p_1, m_width );
    d = straight_seg.get_distance_from_line_seg( &line_seg_ab, near_xy,
        other_near_xy);
    }
else{
    d = get_distance_from_line_seg( &line_seg_ab, near_xy,
        other_near_xy);
    }
AA_DECR_CALL_DEPTH();
return d;
}

double np02_arc::get_distance_from_circle(const np02_circle *c,
    np02_xy *near_xy, np02_xy *circle_near_xy) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AA_ALWAYS_ASSERT(NULL != c);
AUTO_ASSERT(m_radius >= 0.0);
AUTO_ASSERT(c->get_radius() >= 0.0);
double d_from = 0.0;
if( arc_is_approx_straight_line_seg() ){
    np02_line_seg straight_seg;
    straight_seg.init( m_p_0, m_p_1, m_width );
    d_from = straight_seg.get_distance_from_circle(c, near_xy, circle_near_xy );
    }
else{
    np02_xy c_ctr_near_xy;
    const double d = get_distance_from_xy_hw(c->get_ctr(), 0.0, &c_ctr_near_xy);
    const double hw = m_width / 2.0;
    const double& r = c->get_radius();
    d_from = d - (hw + r);
    if( ( NULL != near_xy ) || ( NULL != circle_near_xy ) ){
        const np02_xy lambda_cp = (c->get_ctr()).get_unit_vector_to(c_ctr_near_xy);
        if( NULL != near_xy ){
            near_xy->set_x( c_ctr_near_xy.get_x() - ( hw * lambda_cp.get_x() ));
            near_xy->set_y( c_ctr_near_xy.get_y() - ( hw * lambda_cp.get_y() ));
            }
        if( NULL != circle_near_xy ){
            circle_near_xy->set_x((c->get_ctr()).get_x() + (r * lambda_cp.get_x()));
            circle_near_xy->set_y((c->get_ctr()).get_y() + (r * lambda_cp.get_y()));
            }
        }
    }

AUTO_ASSERT( (NULL == near_xy) || 
    (get_distance_from_xy( *near_xy ) < 
    ( np02_shape::m_small_ratio * (1.0 + m_radius + m_width) )));
AA_DECR_CALL_DEPTH();
return d_from;
}

double np02_arc::get_distance_from_arc(const np02_arc *a,
    np02_xy *near_xy, np02_xy *arc_near_xy) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
double d_from = 0.0;
if( arc_is_approx_straight_line_seg() ){
    np02_line_seg straight_seg;
    straight_seg.init( m_p_0, m_p_1, m_width );
    d_from = straight_seg.get_distance_from_arc(a, near_xy, arc_near_xy );
    AUTO_ASSERT( (NULL == near_xy) || 
        (get_distance_from_xy( *near_xy ) < 
        ( np02_shape::m_small_ratio * (1.0 + m_radius + m_width) )));
    }
else if ( a->arc_is_approx_straight_line_seg()) {
    np02_line_seg straight_seg_a;
    straight_seg_a.init( a->m_p_0, a->m_p_1, a->m_width );
    d_from = get_distance_from_line_seg(&straight_seg_a, near_xy, arc_near_xy );
    AUTO_ASSERT( (NULL == near_xy) || 
        (get_distance_from_xy( *near_xy ) < 
        ( np02_shape::m_small_ratio * (1.0 + m_radius + m_width) )));
    }
else{
    np02_dist_from_xy_xy best_result;
    best_result.m_answer_quality = NP02_ANSWER_QUALITY_LARGE_ERROR;
    best_result.m_distance_from = np02_shape::m_googol;
    best_result.m_near_xy_defined = false;
    best_result.m_other_near_xy_defined = false;
    best_result.m_near_xy.set_xy(0.0,0.0);
    best_result.m_other_near_xy.set_xy(0.0,0.0);
    np02_dist_from_xy_xy result;
    for( int i = 0; i < 6; ++i ){
        switch( i ){
        default:
        case 0:
            get_distance_from_arc_zone_0(a, &result);
            break;
        case 1:
            get_distance_from_arc_zone_1(a, &result);
            break;
        case 2:
            get_distance_from_arc_far_perp_zone(a, &result);
            break;
        case 3:
            a->get_distance_from_arc_zone_0(this, &result);
            result.swap_xy();
            break;
        case 4:
            a->get_distance_from_arc_zone_1(this, &result);
            result.swap_xy();
            break;
        case 5:
            a->get_distance_from_arc_far_perp_zone(this, &result);
            result.swap_xy();
            break;
            }
        if( result.compare(best_result) < 0 ){
            best_result = result;
            }
        }
    
    np02_dist_from_xy_xy result_1;
    get_distance_from_arc_centerline_intersect(a, &result, &result_1);
    if( result.compare(best_result) < 0 ){
        best_result = result;
        }
    if( result_1.compare(best_result) < 0 ){
        best_result = result_1;
        }
    
    if( NULL != near_xy ){
        *near_xy = best_result.m_near_xy;
        }
    if( NULL != arc_near_xy ){
        *arc_near_xy = best_result.m_other_near_xy;
        }
    d_from = best_result.m_distance_from;
    }

AUTO_ASSERT( (NULL == near_xy) || 
    (get_distance_from_xy( *near_xy ) < 
    ( np02_shape::m_small_ratio * (1.0 + m_radius + m_width) )));
AUTO_ASSERT( (NULL == arc_near_xy) || 
    (a->get_distance_from_xy( *arc_near_xy ) < 
    ( np02_shape::m_small_ratio * (1.0 + m_radius + m_width) )));
AA_DECR_CALL_DEPTH();
return d_from;
}

double np02_arc::get_distance_from_line_seg(
    const np02_line_seg *n, np02_xy *near_xy,
    np02_xy *line_seg_near_xy) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AA_ALWAYS_ASSERT(NULL != n);
AUTO_ASSERT(m_radius >= 0.0);
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_1);
double d_from = 0.0;
if( arc_is_approx_straight_line_seg() ){
    np02_line_seg straight_seg;
    straight_seg.init( m_p_0, m_p_1, m_width );
    d_from = straight_seg.get_distance_from_line_seg(n, near_xy, line_seg_near_xy );
    }
else{
    np02_dist_from_xy_xy best_result;
    best_result.m_answer_quality = NP02_ANSWER_QUALITY_LARGE_ERROR;
    best_result.m_distance_from = np02_shape::m_googol;
    best_result.m_near_xy_defined = false;
    best_result.m_other_near_xy_defined = false;
    best_result.m_near_xy.set_xy(0.0,0.0);
    best_result.m_other_near_xy.set_xy(0.0,0.0);
    np02_dist_from_xy_xy result;
    for( int i = 0; i < 6; ++i ){
        switch( i ){
        default:
        case 0:
            get_distance_from_shape_arc_p01(n, true, &result);
            break;
        case 1:
            get_distance_from_shape_arc_p01(n, false, &result);
            break;
        case 2:
            get_distance_from_shape_far_perp_zone(n, &result);
            break;
        case 3:
            get_distance_from_line_seg_far_perp_zone(n, &result);
            break;
        case 4:
            get_distance_from_line_seg_p01(n, false, &result);
            break;
        case 5:
            get_distance_from_line_seg_p01(n, true, &result);
            break;
            }
        if( result.compare(best_result) < 0 ){
            best_result = result;
            }
        }
    
    np02_dist_from_xy_xy result_1;
    get_distance_from_line_seg_centerline_intersect(n, &result, &result_1);
    
    if( result.compare(best_result) < 0 ){
        best_result = result;
        }
    if( result_1.compare(best_result) < 0 ){
        best_result = result_1;
        }
    
    if( NULL != near_xy ){
        *near_xy = best_result.m_near_xy;
        }
    if( NULL != line_seg_near_xy ){
        *line_seg_near_xy = best_result.m_other_near_xy;
        }
    d_from = best_result.m_distance_from;
    }

AUTO_ASSERT( (NULL == near_xy) || 
    (get_distance_from_xy( *near_xy ) < 
    ( np02_shape::m_small_ratio * (1.0 + m_radius + m_width) )));
AA_DECR_CALL_DEPTH();
return d_from;
}

double np02_arc::get_distance_from_rect(const np02_rect *rect,
    np02_xy *near_xy, np02_xy *rect_near_xy) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AA_ALWAYS_ASSERT(NULL != rect);
AUTO_ASSERT(m_radius >= 0.0);
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_1);
double d_from = 0.0;
if( arc_is_approx_straight_line_seg() ){
    np02_line_seg straight_seg;
    straight_seg.init( m_p_0, m_p_1, m_width );
    d_from = straight_seg.get_distance_from_rect(rect, near_xy,
        rect_near_xy );
    }
else{
    np02_dist_from_xy_xy best_result;
    best_result.m_answer_quality = NP02_ANSWER_QUALITY_LARGE_ERROR;
    best_result.m_distance_from = np02_shape::m_googol;
    best_result.m_near_xy_defined = false;
    best_result.m_other_near_xy_defined = false;
    best_result.m_near_xy.set_xy(0.0,0.0);
    best_result.m_other_near_xy.set_xy(0.0,0.0);
    np02_dist_from_xy_xy result;
    for( int i = 0; i < 3; ++i ){
        switch( i ){
        default:
        case 0:
            get_distance_from_shape_arc_p01(rect, true, &result);
            break;
        case 1:
            get_distance_from_shape_arc_p01(rect, false, &result);
            break;
        case 2:
            get_distance_from_shape_far_perp_zone(rect, &result);
            break;
            }
        if( result.compare(best_result) < 0 ){
            best_result = result;
            }
        }
    np02_rect_pt_idx rect_pt_idx = static_cast<np02_rect_pt_idx>(0);
    for( ; rect_pt_idx < NP02_RECT_PT_IDX_COUNT;
        rect_pt_idx = static_cast<np02_rect_pt_idx>(rect_pt_idx + 1) ){
        get_distance_from_rect_pt(rect, rect_pt_idx, &result);
        if( result.compare(best_result) < 0 ){
            best_result = result;
            }
        }
    
    /* if bounding boxes overlap, check rectangle diagonals and parallel bisectors

      A   A   B
      *---+---*  
      |\  |  /|
      | \ | / |
      |  \|/  |
     A+---+---+B
      |  /|\  |
      | / | \ |
      |/  |  \|
     A*---+---*B
          B
    */
    np02_xy rect_bb_xy_min(0.0,0.0);
    np02_xy rect_bb_xy_max(0.0,0.0);
    rect->get_bb( &rect_bb_xy_min, &rect_bb_xy_max );
    const bool bb_overlap = is_bb_overlap( m_bb_xy_min, m_bb_xy_max,
         rect_bb_xy_min, rect_bb_xy_max );
    if( bb_overlap ){
        np02_xy p_a(0.0,0.0);
        np02_xy p_b(0.0,0.0);
        np02_line_seg seg_ab;
        result.m_answer_quality = NP02_ANSWER_QUALITY_SMALL_ERROR;
        result.m_near_xy_defined = true;
        result.m_other_near_xy_defined = true;
        for( int i = 0; i < 4; ++i ){
            switch(i){
            default:
            case 0:
                p_a = rect->get_p00();
                p_b = rect->get_p11();
                break;
            case 1:
                p_a = rect->get_p01();
                p_b = rect->get_p10();
                break;
            case 2:
                p_a = (rect->get_p00()).get_midpoint_to(rect->get_p01());
                p_b = (rect->get_p10()).get_midpoint_to(rect->get_p11());
                break;
            case 3:
                p_a = (rect->get_p00()).get_midpoint_to(rect->get_p10());
                p_b = (rect->get_p01()).get_midpoint_to(rect->get_p11());
                break;
                }
            seg_ab.init( p_a, p_b, 0.0 );
            np02_xy seg_ab_near_xy(0.0,0.0);
            get_distance_from_line_seg( &seg_ab, &(result.m_near_xy), &seg_ab_near_xy );
            result.m_distance_from = rect->get_distance_from_xy(result.m_near_xy,
                &(result.m_other_near_xy));
            if( result.compare(best_result) < 0 ){
                best_result = result;
                }

            /* calculate from arc centerline */
            get_distance_from_xy_hw(seg_ab_near_xy, 0.0, &result.m_near_xy);
            result.m_distance_from = rect->get_distance_from_xy(result.m_near_xy,
                &(result.m_other_near_xy));
            if( result.compare(best_result) < 0 ){
                best_result = result;
                }
            }
        }
    
    if( NULL != near_xy ){
        *near_xy = best_result.m_near_xy;
        }
    if( NULL != rect_near_xy ){
        *rect_near_xy = best_result.m_other_near_xy;
        }
    d_from = best_result.m_distance_from;
    }

AUTO_ASSERT( (NULL == near_xy) || 
    (get_distance_from_xy( *near_xy ) < 
    ( np02_shape::m_small_ratio * (1.0 + m_radius + m_width) )));
AA_DECR_CALL_DEPTH();
return d_from;
}

double np02_arc::get_distance_from_polygon(const np02_polygon *p,
    np02_xy *near_xy, np02_xy *polygon_near_xy) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AA_ALWAYS_ASSERT(NULL != p);
AUTO_ASSERT(m_radius >= 0.0);
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_1);

double d = 0.0;
/* TODO: implement polygon */
AA_ALWAYS_ASSERT(false);

AA_DECR_CALL_DEPTH();
return d;
}

double np02_arc::get_distance_from_shape(const np02_shape *s,
    np02_xy *near_xy, np02_xy *other_near_xy) const{
double d = 0.0;
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
if(NULL != s){
    if( arc_is_approx_straight_line_seg() ){
        np02_line_seg straight_seg;
        straight_seg.init( m_p_0, m_p_1, m_width );
        d = straight_seg.get_distance_from_shape(s, near_xy, other_near_xy );
        }
    else{
        const np02_circle *c = dynamic_cast<const np02_circle *>(s);
        const np02_arc *a = dynamic_cast<const np02_arc *>(s);
        const np02_line_seg *n = dynamic_cast<const np02_line_seg *>(s);
        const np02_rect *r = dynamic_cast<const np02_rect *>(s);
        const np02_polygon *p = dynamic_cast<const np02_polygon *>(s);
        const np02_spline *i = dynamic_cast<const np02_spline *>(s);
        if(NULL != c){ d = get_distance_from_circle(c,near_xy,other_near_xy); }
        else if(NULL!=a){ d=get_distance_from_arc(a, near_xy, other_near_xy); }
        else if(NULL!=n){d=get_distance_from_line_seg(n,near_xy,other_near_xy);}
        else if(NULL!=r){ d=get_distance_from_rect(r, near_xy, other_near_xy);}
        else if(NULL!=p){d=get_distance_from_polygon(p,near_xy,other_near_xy);}
        else if(NULL!=i){ AA_ALWAYS_ASSERT(false);/* spline not supported */}
        else{ AA_ALWAYS_ASSERT(false); }
        }
    }
AA_DECR_CALL_DEPTH();
return d;
}

void np02_arc::translate_no_loc_grid(const np02_xy& dxy){
if( ( dxy.get_x() != 0.0 ) && ( dxy.get_y() != 0.0 ) ){
    init_params p;
    if( arc_is_approx_straight_line_seg() ){
        init3pt_params init3pt_prm;
        const np02_xy p0( m_p_0.get_x() + dxy.get_x(),
                          m_p_0.get_y() + dxy.get_y() );
        const np02_xy p1( m_p_1.get_x() + dxy.get_x(),
                          m_p_1.get_y() + dxy.get_y() );
        init3pt_prm.m_pa = p0;
        init3pt_prm.m_pb.set_x( p0.get_x() + p1.get_x() );
        init3pt_prm.m_pb.set_y( p0.get_y() + p1.get_y() );
        init3pt_prm.m_pc = p1;
        init3pt_prm.m_width = m_width;
        init3pt_prm.m_max_radius = m_radius;
        init3pt_aux_params aux_prm;
        init3pt_to_init_params( init3pt_prm, &p, &aux_prm ); 
        init_force_p01(p, aux_prm );
        }
    else{
        p.m_ctr.set_x( m_ctr.get_x() + dxy.get_x() );
        p.m_ctr.set_y( m_ctr.get_y() + dxy.get_y() );
        p.m_radius = m_radius;
        p.m_start_angle_deg = m_start_angle_deg;
        p.m_end_angle_deg = m_end_angle_deg;
        p.m_width = m_width;    
        init(p);
        }
    }
}

void np02_arc::rotate_no_loc_grid(const np02_xy& rot_ctr,
    const double& rot_deg){
double rot_d = rot_deg;
if( ( rot_d < -360.0 ) || ( rot_d > 360.0 ) ){
    rot_d = fmod( rot_deg, 360.0 );
    }
if( 0.0 != rot_d ){
    init_params p;
    if( arc_is_approx_straight_line_seg() ){
        init3pt_params init3pt_prm;
        np02_xy p0 = m_p_0;
        np02_xy p1 = m_p_1;
        p0.rotate(rot_ctr,rot_deg);
        p1.rotate(rot_ctr,rot_deg);
        init3pt_prm.m_pa = p0;
        init3pt_prm.m_pb.set_x( p0.get_x() + p1.get_x() );
        init3pt_prm.m_pb.set_y( p0.get_y() + p1.get_y() );
        init3pt_prm.m_pc = p1;
        init3pt_prm.m_width = m_width;
        init3pt_prm.m_max_radius = m_radius;
        init3pt_aux_params aux_prm;
        init3pt_to_init_params( init3pt_prm, &p, &aux_prm ); 
        init_force_p01(p, aux_prm );
        }
    else{
        const np02_xy rot_arm_initial(m_ctr.get_x() - rot_ctr.get_x(),
            m_ctr.get_y() - rot_ctr.get_y());
        if( (rot_arm_initial.get_x() != 0.0) || (rot_arm_initial.get_y() != 0.0) ){
            double cos_rot = 1.0;
            double sin_rot = 0.0;
            cos_sin_rot_deg( rot_d, &cos_rot, &sin_rot );
            const np02_xy rot_arm_final(
                (rot_arm_initial.get_x() * cos_rot) -
                (rot_arm_initial.get_y() * sin_rot),
                (rot_arm_initial.get_x() * sin_rot) +
                (rot_arm_initial.get_y() * cos_rot));
            p.m_ctr.set_x( rot_ctr.get_x() + rot_arm_final.get_x() );
            p.m_ctr.set_y( rot_ctr.get_y() + rot_arm_final.get_y() );
            }
        else{
            p.m_ctr = m_ctr;
            }
        p.m_radius = m_radius;
        p.m_start_angle_deg = m_start_angle_deg + rot_d;
        p.m_end_angle_deg = m_end_angle_deg + rot_d;
        if( p.m_start_angle_deg < -360.0 ){
            p.m_start_angle_deg += 360.0;
            p.m_end_angle_deg += 360.0;
            }
        else if( p.m_end_angle_deg > 360.0 ){
            p.m_start_angle_deg -= 360.0;
            p.m_end_angle_deg -= 360.0;
            }
        p.m_width = m_width;    
        init(p);
        }
    }
}

uint64_t np02_arc::hash( const uint64_t& h_in ) const{
uint64_t h = np02_shape::hash( h_in );
h = m_ctr.hash( h );
#if defined( CF01_SUPPORT )
h = cf01_obj_hash( h, m_radius );
h = cf01_obj_hash( h, m_start_angle_deg );
h = cf01_obj_hash( h, m_end_angle_deg );
h = cf01_obj_hash( h, m_width );
#else
h += reinterpret_cast<cf01_uint64>(m_radius);
h ^= ((h << 17) | (h >> 47));
h += reinterpret_cast<cf01_uint64>(m_start_angle_deg);
h ^= ((h << 27) | (h >> 37));
h += reinterpret_cast<cf01_uint64>(m_end_angle_deg);
h ^= ((h << 37) | (h >> 27));
h += reinterpret_cast<cf01_uint64>(m_width);
h ^= ((h << 47) | (h >> 17));
#endif
h = m_p_0.hash( h );
h = m_p_1.hash( h );
h = m_fwd_0.hash( h );
h = m_fwd_1.hash( h );
#if defined( CF01_SUPPORT )
h = cf01_obj_hash( h, m_fwd_dot_0 );
h = cf01_obj_hash( h, m_fwd_dot_1 );
#else
h += reinterpret_cast<cf01_uint64>(m_fwd_dot_0);
h ^= ((h << 27) | (h >> 37));
h += reinterpret_cast<cf01_uint64>(m_fwd_dot_1);
h ^= ((h << 37) | (h >> 27));
#endif
h = m_bb_xy_min.hash( h );
h = m_bb_xy_max.hash( h );
return h;
}

int np02_arc::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const{
int err_cnt = 0;
err_cnt += np02_shape::verify_data(err_msg, err_msg_capacity, err_msg_pos );

if(NP02_SHAPE_TYPE_ARC != get_shape_type() ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "arc: this=%x  m_shape_type=%i\n", this, get_shape_type()); 
    }

const np02_shp_alloc *shp_alloc = get_shp_alloc();
if((NULL != shp_alloc) && 
    (this != shp_alloc->alloc_get_arc_by_idx(get_shape_idx()))){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "arc: this=%x  != (shp_alloc=%x) circle_by_idx(idx=%i)=%x\n",
        this, shp_alloc, get_shape_idx(),
        shp_alloc->alloc_get_arc_by_idx(get_shape_idx())); 
    }

const double big_radius = m_radius + m_width;
static const double little_ratio_error = np02_shape::m_little_ratio;
static const double small_ratio_error = np02_shape::m_small_ratio;
static const double small_ratio_error_sq = small_ratio_error*small_ratio_error;
const double small_distance_error =
    (big_radius > 1.0) ? (small_ratio_error * big_radius) : small_ratio_error;
const double little_distance_error_sq = (big_radius > 1.0) ?
    (little_ratio_error * big_radius * big_radius) : little_ratio_error;
const double small_distance_error_sq = (big_radius > 1.0) ?
    (small_ratio_error * big_radius * big_radius) : small_ratio_error;

/* m_ctr */

if(m_radius < 0.0){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "arc: this=%x m_radius=%f < 0.0\n", this, m_radius);
    }

if( m_start_angle_deg < -720.0 ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "arc: this=%x m_start_angle_deg=%f < -720.0\n", 
        this, m_start_angle_deg);
    }

if( m_start_angle_deg > 720.0 ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "arc: this=%x m_start_angle_deg=%f > 720.0\n", 
        this, m_start_angle_deg);
    }

if( m_end_angle_deg < -720.0 ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "arc: this=%x m_end_angle_deg=%f < -720.0\n", 
        this, m_end_angle_deg);
    }

if( m_end_angle_deg > 720.0 ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "arc: this=%x m_end_angle_deg=%f > 720.0\n", 
        this, m_end_angle_deg);
    }

if( m_end_angle_deg < m_start_angle_deg ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "arc: this=%x end_angle_deg=%f < start_angle_deg=%f\n", 
        this, m_end_angle_deg, m_start_angle_deg);
    }

if( m_end_angle_deg > ( m_start_angle_deg + 360.0 ) ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "arc: this=%x end_angle_deg=%f > (start_angle_deg=%f)+360\n", 
        this, m_end_angle_deg, m_start_angle_deg);
    }

const double start_angle_rad = m_start_angle_deg / m_degrees_per_radian;
const double cos_start_angle = cos(start_angle_rad);
const double sin_start_angle = sin(start_angle_rad);
const np02_xy p_0_check(m_ctr.get_x() + (m_radius*cos_start_angle),  
                        m_ctr.get_y() + (m_radius*sin_start_angle) );
const double p_0_err_sq = m_p_0.get_dsq_to(p_0_check);
if(p_0_err_sq > little_distance_error_sq) {
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "arc: this=%x p_0=(%f,%f) != (%f,%f)="
        "m_ctr(%f,%f) + radius(%f)<(%f)deg  d_err=%f\n", this, 
        m_p_0.get_x(), m_p_0.get_y(), 
        p_0_check.get_x(), p_0_check.get_y(), 
        m_ctr.get_x(), m_ctr.get_y(),
        m_radius, m_start_angle_deg, sqrt(p_0_err_sq) );
    }

const double end_angle_rad = m_end_angle_deg / m_degrees_per_radian;
const double cos_end_angle = cos(end_angle_rad);
const double sin_end_angle = sin(end_angle_rad);
const np02_xy p_1_check(m_ctr.get_x() + (m_radius*cos_end_angle),  
                        m_ctr.get_y() + (m_radius*sin_end_angle) );
const double p_1_err_sq = m_p_1.get_dsq_to(p_1_check);
if(p_1_err_sq > little_distance_error_sq) {
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "arc: this=%x p_1=(%f,%f) != (%f,%f)="
        "m_ctr(%f,%f) + radius(%f)<(%f)deg  d_err=%f\n", this, 
        m_p_1.get_x(), m_p_1.get_y(), 
        p_1_check.get_x(), p_1_check.get_y(), 
        m_ctr.get_x(), m_ctr.get_y(),
        m_radius, m_end_angle_deg, sqrt(p_0_err_sq) );
    }

const double hw = m_width / 2.0;
if(m_width < 0.0){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "arc: this=%x m_width=%f < 0.0\n", this, m_width);
    }

const np02_xy fwd_0_check(-sin_start_angle, cos_start_angle);
const double fwd_0_err_sq = m_fwd_0.get_dsq_to(fwd_0_check);
if(fwd_0_err_sq > small_ratio_error_sq) {
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "arc: this=%x fwd_0=(%f,%f) != (%f,%f)"
        "  start_angle = %fdeg\n", this, 
        m_fwd_0.get_x(), m_fwd_0.get_y(), 
        fwd_0_check.get_x(), fwd_0_check.get_y(), 
        m_start_angle_deg );
    }

const np02_xy fwd_1_check(-sin_end_angle, cos_end_angle);
const double fwd_1_err_sq = m_fwd_1.get_dsq_to(fwd_1_check);
if(fwd_1_err_sq > small_ratio_error_sq) {
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "arc: this=%x fwd_1=(%f,%f) != (%f,%f)"
        "  end_angle = %fdeg\n", this, 
        m_fwd_1.get_x(), m_fwd_1.get_y(), 
        fwd_1_check.get_x(), fwd_1_check.get_y(), 
        m_end_angle_deg );
    }

/* bounding box */
static const double angle_step_degrees = 5.0;
bool bb_check_done = false;
double bb_check_angle_deg = m_start_angle_deg;
while( !bb_check_done ){
    np02_xy p;
    if( bb_check_angle_deg == m_start_angle_deg ){
        p = m_p_0;
        }
    else if( bb_check_angle_deg >= m_end_angle_deg ){
        p = m_p_1;
        bb_check_angle_deg = m_end_angle_deg;
        }
    else{
        const double bb_check_angle_rad =
            bb_check_angle_deg / m_degrees_per_radian;
        p.set_x( m_ctr.get_x() + ( m_radius * cos(bb_check_angle_rad ) ) ); 
        p.set_y( m_ctr.get_y() + ( m_radius * sin(bb_check_angle_rad ) ) ); 
        }
    if( ((p.get_x() - hw) < ( m_bb_xy_min.get_x() - small_distance_error)) ||
        ((p.get_y() - hw) < ( m_bb_xy_min.get_y() - small_distance_error)) ||
        ((p.get_x() + hw) > ( m_bb_xy_max.get_x() + small_distance_error)) ||
        ((p.get_y() + hw) > ( m_bb_xy_max.get_y() + small_distance_error)) ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "arc: this=%x  p(%f,%f)(+/-w/2=%f) angle = %fdeg is outside "
            "bounding box bb_xy_min(%f,%f)/bb_xy_max(%f,%f)\n", this, 
            p.get_x(), p.get_y(), hw, bb_check_angle_deg,
            m_bb_xy_min.get_x(), m_bb_xy_min.get_y(),
            m_bb_xy_max.get_x(), m_bb_xy_max.get_y() );
        }

    if( ( bb_check_angle_deg >= m_end_angle_deg ) || 
        ( bb_check_angle_deg >= ( m_start_angle_deg + 360.0 ) )){
        bb_check_done = true; 
        }
    else{
        bb_check_angle_deg += angle_step_degrees;
        }
    }

return err_cnt;
}

std::ostream& np02_arc::ostream_output(std::ostream& os) const{
os << "<arc>\n";
os << std::hex;
os << "<this>" << this << "</this>\n";
os << std::dec;
np02_shape::ostream_output(os);
os << "<ctr><x>" << m_ctr.get_x() << "</x>";
os << "<y>" << m_ctr.get_y() << "</y></ctr>\n";
os << "<radius>" << m_radius << "</radius>\n";
os << "<start_angle_deg>" << m_start_angle_deg << "</start_angle_deg>\n";
os << "<end_angle_deg>" << m_end_angle_deg << "</end_angle_deg>\n";
os << "<width>" << m_width << "</width>\n";
os << "<p_0><x>" << m_p_0.get_x() << "</x>";
os << "<y>" << m_p_0.get_y() << "</y></p_0>\n";
os << "<p_1><x>" << m_p_1.get_x() << "</x>";
os << "<y>" << m_p_1.get_y() << "</y></p_1>\n";
os << "<fwd_0><x>" << m_fwd_0.get_x() << "</x>";
os << "<y>" << m_fwd_0.get_y() << "</y></fwd_0>\n";
os << "<fwd_1><x>" << m_fwd_1.get_x() << "</x>";
os << "<y>" << m_fwd_1.get_y() << "</y></fwd_1>\n";
os << "<fwd_dot_0>" << m_fwd_dot_0 << "</fwd_dot_0>\n";
os << "<fwd_dot_1>" << m_fwd_dot_1 << "</fwd_dot_1>\n";
os << "<bb_xy_min><x>" << m_bb_xy_min.get_x() << "</x>";
os << "<y>" << m_bb_xy_min.get_y() << "</y></bb_xy_min>\n";
os << "<bb_xy_max><x>" << m_bb_xy_max.get_x() << "</x>";
os << "<y>" << m_bb_xy_max.get_y() << "</y></bb_xy_max>\n";
os << "</arc>\n";
return os;
}

void np02_arc::write_bmp_file(const np02_xy& xy_min,
    const double& pixel_num, const np02_bmp_color& color,
    np02_bmp_file *bmp_file) const{
np02_shape::write_bmp_file(xy_min, pixel_num, color, bmp_file);
}

void np02_arc::write_dxf_file(const std::string& layer,
    const uint8_t& color, np02_dxf_file *dxf_file) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AUTO_ASSERT(NULL != dxf_file);
/* centerline */
dxf_file->draw_arc( layer, m_ctr.get_x(), m_ctr.get_y(), m_radius,
    m_start_angle_deg, m_end_angle_deg, color );

/* boundary */
if( m_width > 0.0 ){
    const double hw = m_width / 2.0;
    const double big_r = m_radius + hw;
    dxf_file->draw_arc( layer, m_ctr.get_x(), m_ctr.get_y(), big_r,
        m_start_angle_deg, m_end_angle_deg, color );
    if( m_radius > hw ){
        const double small_r = m_radius - hw;
        dxf_file->draw_arc( layer, m_ctr.get_x(), m_ctr.get_y(), small_r,
            m_start_angle_deg, m_end_angle_deg, color );
        }
    dxf_file->draw_arc( layer, m_p_0.get_x(), m_p_0.get_y(), hw,
        m_start_angle_deg-180.0, m_start_angle_deg, color );
    dxf_file->draw_arc( layer, m_p_1.get_x(), m_p_1.get_y(), hw,
        m_end_angle_deg, m_end_angle_deg+180.0, color );
    }

AA_DECR_CALL_DEPTH();
}

size_t np02_arc::arc_circle_centerline_intersect(const np02_circle *c,
    np02_xy *xy_intsct_0, np02_xy *xy_intsct_1) const{
size_t intsct_count = 0;
if( NULL == c ){
    intsct_count = 0;
    }
else if( ( m_ctr == c->get_ctr() ) && ( m_radius == c->get_radius() ) ){
    if( m_p_0 == m_p_1 ){
        intsct_count = 1;
        if( NULL != xy_intsct_0 ){ *xy_intsct_0 = m_p_0; }
        }
    else{
        intsct_count = std::numeric_limits<size_t>::max();
        if( NULL != xy_intsct_0 ){ *xy_intsct_0 = m_p_0; }
        if( NULL != xy_intsct_1 ){ *xy_intsct_1 = m_p_1; }
        }
    }
else{
    np02_xy ccxy0(0.0,0.0);
    np02_xy ccxy1(0.0,0.0);
    const size_t circle_circle_intsct_count = 
        np02_circle::circle_circle_intsct( m_ctr, m_radius, c->get_ctr(),
            c->get_radius(), &ccxy0, &ccxy1 );
    if( ( circle_circle_intsct_count > 0 ) && ( is_in_perp_zone(ccxy0) ) ){
        ++intsct_count;
        if( NULL != xy_intsct_0 ){ *xy_intsct_0 = ccxy0; }        
        }
    if( ( circle_circle_intsct_count > 0 ) && ( is_in_perp_zone(ccxy1) ) ){
        ++intsct_count;
        if( NULL != xy_intsct_1 ){ *xy_intsct_1 = ccxy1; }        
        }
    }

return intsct_count;
}

int np02_arc::verify_data_arc_circle_centerline_intersect_result(
    const np02_circle *c, np02_xy *xy_intsct_0, np02_xy *xy_intsct_1,
    const size_t& intsct_count, char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;

const double& r_a = m_radius;
const double& r_b = c->get_radius();
const np02_xy& ctr_a = m_ctr;
const np02_xy& ctr_b = c->get_ctr();
const double r_a_sq = r_a * r_a;
const double r_b_sq = r_b * r_b;
const double r_ab_sq_sum = r_a_sq + r_b_sq;
const double r_ab = fabs(r_a) + fabs(r_b);
const double max_d_err = ( r_ab > 1.0 ) ? 
    ( r_ab * m_small_ratio ) : m_small_ratio;
const double max_d_err_sq = ( r_ab_sq_sum > 1.0 ) ? 
    ( r_ab_sq_sum * m_little_ratio_sq ) : m_little_ratio_sq;

if( ( intsct_count > 0 ) && ( NULL != xy_intsct_0 ) ){
    const double d_in_perp_a_0 = distance_in_perp_zone(*xy_intsct_0);
    if( d_in_perp_a_0 < -max_d_err ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_in_perp_a_0=%g < 0.0\n", d_in_perp_a_0 ); 
        }
 
    const double dsq_a_0 = xy_intsct_0->get_dsq_to( ctr_a );
    const double dsq_a_0_err = fabs( dsq_a_0 - r_a_sq );
    if( dsq_a_0_err > max_d_err_sq ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_a_0=%g != r_a=%g\n", sqrt(dsq_a_0), r_a ); 
        }

    const double dsq_b_0 = xy_intsct_0->get_dsq_to( ctr_b );
    const double dsq_b_0_err = fabs( dsq_b_0 - r_b_sq );
    if( dsq_b_0_err > max_d_err_sq ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_b_0=%g != r_b=%g\n", sqrt(dsq_b_0), r_b ); 
        }
    }

if( ( intsct_count > 1 ) && ( NULL != xy_intsct_1 ) ){
    const double d_in_perp_a_1 = distance_in_perp_zone(*xy_intsct_1);
    if( d_in_perp_a_1 < -max_d_err ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_in_perp_a_1=%g < 0.0\n", d_in_perp_a_1 ); 
        }

    const double dsq_a_1 = xy_intsct_1->get_dsq_to( ctr_a );
    const double dsq_a_1_err = fabs( dsq_a_1 - r_a_sq );
    if( dsq_a_1_err > max_d_err_sq ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_a_1=%g != r_a=%g\n", sqrt(dsq_a_1), r_a ); 
        }

    const double dsq_b_1 = xy_intsct_1->get_dsq_to( ctr_b );
    const double dsq_b_1_err = fabs( dsq_b_1 - r_b_sq );
    if( dsq_b_1_err > max_d_err_sq ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_b_1=%g != r_b=%g\n", sqrt(dsq_b_1), r_b ); 
        }
    }

if( err_cnt > 0 ){    
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "verify_data_arc_circle_centerline_intersect_result( ctr_a=[%g,%g],"
        " r_a=%g, ctr_b=[%g,%g], r_b=%g, xy_intsct_0=[%g,%g],"
        " xy_intsct_1=[%g,%g], intsct_count=%i\n",
        ctr_a.get_x(), ctr_a.get_y(), r_a, 
        ctr_b.get_x(), ctr_b.get_y(), r_b,
        ( NULL == xy_intsct_0 ) ? 0.0 : xy_intsct_0->get_x(),
        ( NULL == xy_intsct_0 ) ? 0.0 : xy_intsct_0->get_y(),
        ( NULL == xy_intsct_1 ) ? 0.0 : xy_intsct_1->get_x(),
        ( NULL == xy_intsct_1 ) ? 0.0 : xy_intsct_1->get_y(),
        intsct_count );
    }

return err_cnt;
}

size_t np02_arc::arc_arc_centerline_intersect(const np02_arc *a,
    np02_xy *xy_intsct_0, np02_xy *xy_intsct_1) const{
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
size_t intsct_count = 0;
if( NULL == a ){
    intsct_count = 0;
    }
else if( arc_is_approx_straight_line_seg() ){
    np02_line_seg straight_seg;
    straight_seg.init( m_p_0, m_p_1, 0.0 );
    intsct_count = a->arc_seg_centerline_intersect( &straight_seg,
        xy_intsct_0, xy_intsct_1);
    }
else if( a->arc_is_approx_straight_line_seg() ){
    np02_line_seg straight_seg_a;
    straight_seg_a.init( a->get_p_0(), a->get_p_1(), 0.0 );
    intsct_count = arc_seg_centerline_intersect( &straight_seg_a,
        xy_intsct_0, xy_intsct_1);
    }
else if( ( m_ctr == a->get_ctr() ) && ( m_radius == a->get_radius() ) ){
    if( 0.0 == m_radius ){
        intsct_count = 1;
        if( NULL != xy_intsct_0 ){ *xy_intsct_0 = m_ctr; }
        }
    else{
        const bool this_z_a0 = is_in_perp_zone(a->m_p_0);
        const bool this_z_a1 = is_in_perp_zone(a->m_p_1);
        const bool a_z_0 = a->is_in_perp_zone(m_p_0);
        const bool a_z_1 = a->is_in_perp_zone(m_p_1);
        if( (!this_z_a0) && (!this_z_a1) && (!a_z_0) && (!a_z_1) ){
            /*     0t....
                  tt     .a1
                t           a
              1b      +      a
                :           a
                  ..     .a0
                    .....
            */
            intsct_count = 0;
            }
        else if( this_z_a0 && this_z_a1 && a_z_0 && a_z_1 ){
            if( m_p_0 == m_p_1 ){
                 if ( a->m_p_0 == a->m_p_1 ) {
                    if ( m_p_0 == a->m_p_0 ) {
                        /*      .....
                              ..   01ta01
                            :           :
                           :      +      :
                            :           :
                              ..     ..
                                .....
                        */
                        intsct_count = 1;
                        if( NULL != xy_intsct_0 ){ *xy_intsct_0 = m_p_0; }
                        }
                    else{
                        /* unexpected */
                        intsct_count = 1;
                        if( NULL != xy_intsct_0 ){ *xy_intsct_0 = a->m_p_0; }
                        }
                    }
                else{
                    /* unexpected */
                    intsct_count = 1;
                    if( NULL != xy_intsct_0 ){ *xy_intsct_0 = m_p_0; }
                    }
                }
            else if ( a->m_p_0 == a->m_p_1 ) {
                /* unexpected */
                intsct_count = 1;
                if( NULL != xy_intsct_0 ){ *xy_intsct_0 = a->m_p_0; }
                }
            else if ( ( m_p_0 == a->m_p_1 ) && ( m_p_1 == a->m_p_0 ) ){
                /*      ttttt
                      tt    0ta1
                    t           a
                   t      +      a
                    t            a
                      tt    1ta0
                        ttttt
                */
                intsct_count = 2;
                if( NULL != xy_intsct_0 ){ *xy_intsct_0 = m_p_0; }
                if( NULL != xy_intsct_1 ){ *xy_intsct_1 = m_p_1; }
                }
            else{
                /*      ttttt
                      tt     ta1
                    t         0ta
                   t      +      a
                    t          1ta
                      tt     ta0
                        ttttt
                */
                intsct_count = std::numeric_limits<size_t>::max();
                if( NULL != xy_intsct_0 ){ *xy_intsct_0 = m_p_0; }
                if( NULL != xy_intsct_1 ){ *xy_intsct_1 = a->m_p_0; }
                }
            }
        else if( this_z_a0 && (!this_z_a1) && (!a_z_0) && a_z_1 ){
            /*     0t....
                  tt     .a1
                t           a
               t      +    1ta
                t          ta
                  tt     ta0
                    ttttt
            */
            if( NULL != xy_intsct_0 ){
                *xy_intsct_0 = m_p_1;
                }
            if( m_p_1 == a->m_p_0 ){
                intsct_count = 1;
                }
            else{
                if( NULL != xy_intsct_1 ){ *xy_intsct_1 = a->m_p_0; }
                intsct_count = std::numeric_limits<size_t>::max();
                }
            }
        else if( (!this_z_a0) && this_z_a1 && a_z_0 && (!a_z_1) ){
            /*      ttttt
                  tt     ta1
                t          ta
              1t      +    0ta
                :           a
                  ..     .a0
                    .....
            */
            if( NULL != xy_intsct_0 ){
                *xy_intsct_0 = m_p_0;
                }
            if( m_p_0 == a->m_p_1 ){
                intsct_count = 1;
                }
            else{
                if( NULL != xy_intsct_1 ){ *xy_intsct_1 = a->m_p_1; }
                intsct_count = std::numeric_limits<size_t>::max();
                }
            }
        else if( m_p_0 == m_p_1){
            if( a_z_0 ){
                /*      .....
                      ..     .a1
                    :           a
                   :      +   01ta
                    :           a
                      ..     .a0
                        .....
                */
                intsct_count = 1;
                if( NULL != xy_intsct_0 ){ *xy_intsct_0 = m_p_0; }
                }
            else{
                /*      .....
                      ..     .a1
                    :           a
                 01t      +      a
                    :           a
                      ..     .a0
                        .....
                */
                intsct_count = 0;
                }
            }
        else if( a->m_p_0 == a->m_p_1){
            if( this_z_a0 ){
                /*      .....
                      ..     1t
                    :           t
                   :      +     ta01
                    :           t
                      ..     0t
                        .....
                */
                intsct_count = 1;
                if( NULL != xy_intsct_0 ){ *xy_intsct_0 = a->m_p_0; }
                }
            else{
                /*      .....
                      t.     ..
                    t           :
                 01t      +      a01
                    t           :
                      t.     ..
                        .....
                */
                intsct_count = 0;
                }
            }
        else{
            intsct_count = 0; /* unexpected */
            }      
        }
    }
else{
    np02_xy ccxy0(0.0,0.0);
    np02_xy ccxy1(0.0,0.0);
    const size_t circle_circle_intsct_count = 
        np02_circle::circle_circle_intsct( m_ctr, m_radius, a->get_ctr(),
            a->get_radius(), &ccxy0, &ccxy1 );
    if( ( circle_circle_intsct_count > 0 ) && ( is_in_perp_zone(ccxy0) ) &&
        ( a->is_in_perp_zone(ccxy0) ) ){
        ++intsct_count;
        if( NULL != xy_intsct_0 ){ *xy_intsct_0 = ccxy0; }        
        }
    if( ( circle_circle_intsct_count > 1 ) && ( is_in_perp_zone(ccxy1) ) &&
        ( a->is_in_perp_zone(ccxy1) ) ){
        ++intsct_count;
        if( NULL != xy_intsct_1 ){
            if( 1 ==  intsct_count ){
                *xy_intsct_0 = ccxy1;
                }
            else{
                *xy_intsct_1 = ccxy1;
                }
            }        
        }
    }

AUTO_ASSERT( 0 == verify_data_arc_arc_centerline_intersect_result( a,
    xy_intsct_0, xy_intsct_1, intsct_count,
    AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ) );
return intsct_count;
}

int np02_arc::verify_data_arc_arc_centerline_intersect_result( const np02_arc *b,
    np02_xy *xy_intsct_0, np02_xy *xy_intsct_1, const size_t& intsct_count,
    char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const{
int err_cnt = 0;

const double& r_a = m_radius;
const double& r_b = b->get_radius();
const np02_xy& ctr_a = m_ctr;
const np02_xy& ctr_b = b->get_ctr();
const double r_a_sq = r_a * r_a;
const double r_b_sq = r_b * r_b;
const double r_ab_sq_sum = r_a_sq + r_b_sq;
const double r_ab = fabs(r_a) + fabs(r_b);
const double max_d_err = ( r_ab > 1.0 ) ? 
    ( r_ab * m_small_ratio ) : m_small_ratio;
const double max_d_err_sq = ( r_ab_sq_sum > 1.0 ) ? 
    ( r_ab_sq_sum * m_little_ratio_sq ) : m_little_ratio_sq;

if( ( intsct_count > 0 ) && ( NULL != xy_intsct_0 ) ){
    const double d_in_perp_a_0 = distance_in_perp_zone(*xy_intsct_0);
    if( d_in_perp_a_0 < -max_d_err ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_in_perp_a_0=%g < 0.0\n", d_in_perp_a_0 ); 
        }
 
    const double d_in_perp_b_0 = b->distance_in_perp_zone(*xy_intsct_0);
    if( d_in_perp_b_0 < -max_d_err ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_in_perp_b_0=%g < 0.0\n", d_in_perp_b_0 ); 
        }
 
    const double dsq_a_0 = xy_intsct_0->get_dsq_to( ctr_a );
    const double dsq_a_0_err = fabs( dsq_a_0 - r_a_sq );
    if( dsq_a_0_err > max_d_err_sq ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_a_0=%g != r_a=%g\n", sqrt(dsq_a_0), r_a ); 
        }

    const double dsq_b_0 = xy_intsct_0->get_dsq_to( ctr_b );
    const double dsq_b_0_err = fabs( dsq_b_0 - r_b_sq );
    if( dsq_b_0_err > max_d_err_sq ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_b_0=%g != r_b=%g\n", sqrt(dsq_b_0), r_b ); 
        }
    }

if( ( intsct_count > 1 ) && ( NULL != xy_intsct_1 ) ){
    const double d_in_perp_a_1 = distance_in_perp_zone(*xy_intsct_1);
    if( d_in_perp_a_1 < -max_d_err ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_in_perp_a_1=%g < 0.0\n", d_in_perp_a_1 ); 
        }

    const double d_in_perp_b_1 = distance_in_perp_zone(*xy_intsct_1);
    if( d_in_perp_b_1 < -max_d_err ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_in_perp_b_1=%g < 0.0\n", d_in_perp_b_1 ); 
        }

    const double dsq_a_1 = xy_intsct_1->get_dsq_to( ctr_a );
    const double dsq_a_1_err = fabs( dsq_a_1 - r_a_sq );
    if( dsq_a_1_err > max_d_err_sq ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_a_1=%g != r_a=%g\n", sqrt(dsq_a_1), r_a ); 
        }

    const double dsq_b_1 = xy_intsct_1->get_dsq_to( ctr_b );
    const double dsq_b_1_err = fabs( dsq_b_1 - r_b_sq );
    if( dsq_b_1_err > max_d_err_sq ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_b_1=%g != r_b=%g\n", sqrt(dsq_b_1), r_b ); 
        }
    }

if( err_cnt > 0 ){    
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "verify_data_arc_arc_centerline_intersect_result( ctr_a=[%g,%g],"
        " r_a=%g, ctr_b=[%g,%g], r_b=%g, xy_intsct_0=[%g,%g],"
        " xy_intsct_1=[%g,%g], intsct_count=%i\n",
        ctr_a.get_x(), ctr_a.get_y(), r_a, 
        ctr_b.get_x(), ctr_b.get_y(), r_b,
        ( NULL == xy_intsct_0 ) ? 0.0 : xy_intsct_0->get_x(),
        ( NULL == xy_intsct_0 ) ? 0.0 : xy_intsct_0->get_y(),
        ( NULL == xy_intsct_1 ) ? 0.0 : xy_intsct_1->get_x(),
        ( NULL == xy_intsct_1 ) ? 0.0 : xy_intsct_1->get_y(),
        intsct_count );
    }

return err_cnt;
}

size_t np02_arc::arc_seg_centerline_intersect(const np02_line_seg *n,
    np02_xy *xy_intsct_0, np02_xy *xy_intsct_1) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AA_ALWAYS_ASSERT(NULL != n);
AUTO_ASSERT(m_radius >= 0.0);
size_t intsct_count = 0;

if( arc_is_approx_straight_line_seg() ){
    np02_line_seg straight_seg;
    straight_seg.init( m_p_0, m_p_1, 0.0 );
    intsct_count = straight_seg.seg_seg_centerline_intersect( n,
        xy_intsct_0, xy_intsct_1);
    }
else{
    /* find near point of center on infinite line */
    const np02_xy& line_seg_fwd = n->get_fwd();
    const double fwd_cross_ctr = line_seg_fwd.cross(m_ctr);
    const double d = fabs( n->get_fwd_cross_01() - fwd_cross_ctr );
    if( d <= m_radius ){
        /* find intersections of circle with infinite line */
        const double& fwd_x = line_seg_fwd.get_x();
        const double& fwd_y = line_seg_fwd.get_y();
        const double& fwd_cross_01 = n->get_fwd_cross_01();
        const double fwd_dot_ctr = line_seg_fwd.dot(m_ctr);
        const double d_sq = d * d;
        const double r_sq = m_radius * m_radius;
        const double h_sq = r_sq - d_sq;
        AUTO_ASSERT(h_sq >= 0.0);
        const double h = ( h_sq > 0.0 ) ? sqrt( h_sq ) : 0.0;
        const int i_end = (h > 0.0 ) ? 2 : 1;
        for( int i = 0; i < i_end; ++i ){
            const double fwd_dot_p = ( 0 == i ) ?
                ( fwd_dot_ctr - h ) : ( fwd_dot_ctr + h );
            if( ( n->get_fwd_dot_0() <= fwd_dot_p ) &&
                ( fwd_dot_p <= n->get_fwd_dot_1() ) ){
                /* intersection P lies along line segment */
                const np02_xy p(( fwd_x * fwd_dot_p ) - ( fwd_y * fwd_cross_01 ),
                                ( fwd_y * fwd_dot_p ) + ( fwd_x * fwd_cross_01 ));
                if( is_in_perp_zone(p) ){
                    /* intersection P lies along arc */
                    ++intsct_count;
                    if( 1 == intsct_count ){
                        if( NULL != xy_intsct_0 ){ *xy_intsct_0 = p; }
                        }
                    else{
                        if( NULL != xy_intsct_1 ){ *xy_intsct_1 = p; }
                        }
                    }
                }
            }
        }
    }
AA_DECR_CALL_DEPTH();
AUTO_ASSERT( 0 == verify_data_arc_seg_centerline_intersect_result( n,
    xy_intsct_0, xy_intsct_1, intsct_count,
    AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ) );
return intsct_count;
}

int np02_arc::verify_data_arc_seg_centerline_intersect_result(
    const np02_line_seg *n, np02_xy *xy_intsct_0, np02_xy *xy_intsct_1,
    const size_t& intsct_count, char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;

const double& r_a = m_radius;
const np02_xy& ctr_a = m_ctr;
const double seg_len = (n->get_fwd_dot_1()) - (n->get_fwd_dot_0());
const double seg_len_sq = seg_len*seg_len;
const double r_a_sq = r_a * r_a;
const double len_sq_sum = r_a_sq + seg_len_sq;
const double len_ab = fabs(r_a) + seg_len;
const double max_d_err = ( len_ab > 1.0 ) ? 
    ( len_ab * m_small_ratio ) : m_small_ratio;
const double max_d_err_sq = ( len_sq_sum > 1.0 ) ? 
    ( len_sq_sum * m_little_ratio_sq ) : m_little_ratio_sq;

if( ( intsct_count > 0 ) && ( NULL != xy_intsct_0 ) ){
    const double d_in_perp_a_0 = distance_in_perp_zone(*xy_intsct_0);
    if( d_in_perp_a_0 < -max_d_err ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_in_perp_a_0=%g < 0.0\n", d_in_perp_a_0 ); 
        }
 
    const double dsq_a_0 = xy_intsct_0->get_dsq_to( ctr_a );
    const double dsq_a_0_err = fabs( dsq_a_0 - r_a_sq );
    if( dsq_a_0_err > max_d_err_sq ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_a_0=%g != r_a=%g\n", sqrt(dsq_a_0), r_a ); 
        }

    const double d_from_seg_0 = n->get_distance_from_xy_hw(
        *xy_intsct_0, 0.0, NULL );
    if( fabs(d_from_seg_0) > max_d_err ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_from_seg_0=%g != 0.0\n", d_from_seg_0 ); 
        }
    }

if( ( intsct_count > 1 ) && ( NULL != xy_intsct_1 ) ){
    const double d_in_perp_a_1 = distance_in_perp_zone(*xy_intsct_1);
    if( d_in_perp_a_1 < -max_d_err ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_in_perp_a_1=%g < 0.0\n", d_in_perp_a_1 ); 
        }
 
    const double dsq_a_1 = xy_intsct_1->get_dsq_to( ctr_a );
    const double dsq_a_1_err = fabs( dsq_a_1 - r_a_sq );
    if( dsq_a_1_err > max_d_err_sq ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_a_1=%g != r_a=%g\n", sqrt(dsq_a_1), r_a ); 
        }

    const double d_from_seg_1 = n->get_distance_from_xy_hw(
        *xy_intsct_1, 0.0, NULL );
    if( fabs(d_from_seg_1) > max_d_err ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d_from_seg_1=%g != 0.0\n", d_from_seg_1 ); 
        }
    }

if( err_cnt > 0 ){
    const np02_xy& seg_p0 = n->get_p_0();
    const np02_xy& seg_p1 = n->get_p_1();
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "verify_data_arc_seg_centerline_intersect_result( ctr_a=[%g,%g],"
        " r_a=%g, seg_p0=[%g,%g], seg_p1=[%g,%g], xy_intsct_0=[%g,%g],"
        " xy_intsct_1=[%g,%g], intsct_count=%i\n",
        ctr_a.get_x(), ctr_a.get_y(), r_a, 
        seg_p0.get_x(), seg_p0.get_y(), seg_p1.get_x(), seg_p1.get_y(),
        ( NULL == xy_intsct_0 ) ? 0.0 : xy_intsct_0->get_x(),
        ( NULL == xy_intsct_0 ) ? 0.0 : xy_intsct_0->get_y(),
        ( NULL == xy_intsct_1 ) ? 0.0 : xy_intsct_1->get_x(),
        ( NULL == xy_intsct_1 ) ? 0.0 : xy_intsct_1->get_y(),
        intsct_count );
    }

return err_cnt;
}

void np02_arc::init_bb(){
const double hw = m_width / 2.0;
const double big_radius = m_radius + hw;
if( m_p_0.get_x() < m_p_1.get_x() ) {
    m_bb_xy_min.set_x( m_p_0.get_x() - hw );
    m_bb_xy_max.set_x( m_p_1.get_x() + hw );
    }
else{
    m_bb_xy_min.set_x( m_p_1.get_x() - hw );
    m_bb_xy_max.set_x( m_p_0.get_x() + hw );
    }
if( m_p_0.get_y() < m_p_1.get_y() ) {
    m_bb_xy_min.set_y( m_p_0.get_y() - hw );
    m_bb_xy_max.set_y( m_p_1.get_y() + hw );
    }
else{
    m_bb_xy_min.set_y( m_p_1.get_y() - hw );
    m_bb_xy_max.set_y( m_p_0.get_y() + hw );
    }
if( ( ( m_start_angle_deg <=    0.0 ) && (    0.0 <= m_end_angle_deg ) ) ||
    ( ( m_start_angle_deg <=  360.0 ) && (  360.0 <= m_end_angle_deg ) ) ||
    ( ( m_start_angle_deg <=  720.0 ) && (  720.0 <= m_end_angle_deg ) ) ){
    const double bb_x_max = m_ctr.get_x() + big_radius;
    if( bb_x_max > m_bb_xy_max.get_x() ) { m_bb_xy_max.set_x(bb_x_max); }
    }
if( ( ( m_start_angle_deg <= -270.0 ) && ( -270.0 <= m_end_angle_deg ) ) ||
    ( ( m_start_angle_deg <=   90.0 ) && (   90.0 <= m_end_angle_deg ) ) ||
    ( ( m_start_angle_deg <=  450.0 ) && (  450.0 <= m_end_angle_deg ) ) ){
    const double bb_y_max = m_ctr.get_y() + big_radius;
    if( bb_y_max > m_bb_xy_max.get_y() ) { m_bb_xy_max.set_y(bb_y_max); }
    }
if( ( ( m_start_angle_deg <= -180.0 ) && ( -180.0 <= m_end_angle_deg ) ) ||
    ( ( m_start_angle_deg <=  180.0 ) && (  180.0 <= m_end_angle_deg ) ) ||
    ( ( m_start_angle_deg <=  540.0 ) && (  540.0 <= m_end_angle_deg ) ) ){
    const double bb_x_min = m_ctr.get_x() - big_radius;
    if( bb_x_min < m_bb_xy_min.get_x() ) { m_bb_xy_min.set_x(bb_x_min); }
    }
if( ( ( m_start_angle_deg <=  -90.0 ) && (  -90.0 <= m_end_angle_deg ) ) ||
    ( ( m_start_angle_deg <=  270.0 ) && (  270.0 <= m_end_angle_deg ) ) ||
    ( ( m_start_angle_deg <=  630.0 ) && (  630.0 <= m_end_angle_deg ) ) ){
    const double bb_y_min = m_ctr.get_y() - big_radius;
    if( bb_y_min < m_bb_xy_min.get_y() ) { m_bb_xy_min.set_y(bb_y_min); }
    }
}

double np02_arc::get_distance_from_xy_hw(const np02_xy& xy,
        const double& hw, np02_xy *near_xy) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AUTO_ASSERT(m_radius >= 0.0);
AUTO_ASSERT(hw >= 0.0);
double d = 0.0;

if( arc_is_approx_straight_line_seg() ){
    np02_line_seg straight_seg;
    straight_seg.init( m_p_0, m_p_1, m_width );
    d = straight_seg.get_distance_from_xy_hw(xy, hw, near_xy );
    }
else{
    bool in_perp_zone = is_in_perp_zone(xy);
    if( in_perp_zone ){
        const double ctr_dx = xy.get_x() - m_ctr.get_x();
        const double ctr_dy = xy.get_y() - m_ctr.get_y();
        const double ctr_dsq = (ctr_dx*ctr_dx) + (ctr_dy*ctr_dy);
        const double d_to_ctr = sqrt(ctr_dsq);
        const double d_arc_ctr = d_to_ctr - m_radius;
    
        /*
                           outside                 .
                                                 .
                                 ***           .
          .                 *           *    .
              .          *                 +P0
                  .    *      inside      .      
                      +P1               .
                           .           .        
                                .    .
                                   +ctr
        */
        if( d_arc_ctr > 0.0 ){
            /* outside */
            d = d_arc_ctr - hw;
            if(NULL != near_xy){
                const double radius_big = m_radius + hw;
                if(radius_big < np02_shape::m_teensy_ratio){
                    *near_xy = m_ctr; 
                    }
                else if(d_to_ctr < np02_shape::m_teensy_ratio){
                    near_xy->set_x(m_ctr.get_x()+(radius_big * m_fwd_0.get_y()));
                    near_xy->set_y(m_ctr.get_y()-(radius_big * m_fwd_0.get_x()));
                    }
                else{ 
                    const double f_out = radius_big / d_to_ctr;
                    near_xy->set_x(m_ctr.get_x() + (f_out * ctr_dx));
                    near_xy->set_y(m_ctr.get_y() + (f_out * ctr_dy));
                    }
                AUTO_ASSERT( 
                    fabs(fabs( near_xy->get_distance_to( m_ctr ) - m_radius ) - hw)
                    <= ( np02_shape::m_small_ratio * (1.0 + radius_big) ) );
                }
            }
        else{
            /* inside */
            d = -d_arc_ctr - hw;
            np02_xy nrxy(0.0,0.0);
            const double radius_small = (m_radius > hw) ? (m_radius-hw) : 0.0;
            if(radius_small < np02_shape::m_teensy_ratio){
                nrxy = m_ctr;
                AUTO_ASSERT( in_perp_zone );
                }
            else{
                if(d_to_ctr < np02_shape::m_teensy_ratio){
                    nrxy.set_x(m_ctr.get_x()+(radius_small * m_fwd_0.get_y()));
                    nrxy.set_y(m_ctr.get_y()-(radius_small * m_fwd_0.get_x()));
                    }
                else{ 
                    const double f_in = radius_small / d_to_ctr;
                    nrxy.set_x(m_ctr.get_x() + (f_in * ctr_dx));
                    nrxy.set_y(m_ctr.get_y() + (f_in * ctr_dy));
                    }
                in_perp_zone = is_in_perp_zone(nrxy);
                }

            if( in_perp_zone && ( NULL != near_xy ) ){
                *near_xy = nrxy;
                AUTO_ASSERT( 
                    fabs(fabs( near_xy->get_distance_to( m_ctr ) - m_radius ) - hw)
                    <= ( np02_shape::m_small_ratio * (1.0 + m_radius + hw) ) );
                }
            }
        }
    if( !in_perp_zone ){
        const double dx_p0 = xy.get_x() - m_p_0.get_x();
        const double dy_p0 = xy.get_y() - m_p_0.get_y();
        const double dsq_p0 = (dx_p0*dx_p0) + (dy_p0*dy_p0);
        const double dx_p1 = xy.get_x() - m_p_1.get_x();
        const double dy_p1 = xy.get_y() - m_p_1.get_y();
        const double dsq_p1 = (dx_p1*dx_p1) + (dy_p1*dy_p1);
        if( dsq_p0 <= dsq_p1 ){
            /* p0 closest */
            const double d_to_p0 = sqrt(dsq_p0);
            d = d_to_p0 - hw;
            if(NULL != near_xy){
                if(hw < np02_shape::m_teensy_ratio) {
                    *near_xy = m_p_0; 
                }
                else if(d_to_p0 < np02_shape::m_teensy_ratio){
                    near_xy->set_x(m_p_0.get_x() - (hw * m_fwd_0.get_x()));
                    near_xy->set_y(m_p_0.get_y() - (hw * m_fwd_0.get_y()));
                }
                else{
                    const double f0 = hw / d_to_p0;
                    near_xy->set_x(m_p_0.get_x() + (f0 * dx_p0));
                    near_xy->set_y(m_p_0.get_y() + (f0 * dy_p0));
                    }
    
                AUTO_ASSERT( fabs( near_xy->get_distance_to( m_p_0 ) - hw)
                    <= ( np02_shape::m_small_ratio * (1.0 + m_radius + hw ) ) );
                }
            }
        else{
            /* p1 closest */
            const double d_to_p1 = sqrt(dsq_p1);
            d = d_to_p1 - hw;
            if(NULL != near_xy){
                if(hw < np02_shape::m_teensy_ratio) {
                    *near_xy = m_p_1; 
                    }
                else if(d_to_p1 < np02_shape::m_teensy_ratio){
                    near_xy->set_x(m_p_1.get_x() + (hw * m_fwd_1.get_x()));
                    near_xy->set_y(m_p_1.get_y() + (hw * m_fwd_1.get_y()));
                    }
                else{
                    const double f1 = hw / d_to_p1;
                    near_xy->set_x(m_p_1.get_x() + (f1 * dx_p1));
                    near_xy->set_y(m_p_1.get_y() + (f1 * dy_p1));
                    }
                AUTO_ASSERT( fabs( near_xy->get_distance_to( m_p_1 ) - hw)
                    <= ( np02_shape::m_small_ratio * (1.0 + m_radius + hw ) ) );
                }
            }
        }
    }

AA_DECR_CALL_DEPTH();
return d;
}

/*                                      
                                                   P
                                                  ^
                     ****         zone 0         ^
                    **C***                      ^   perp
                     ***.***                   ^    zone       ****
                      ***.***              ** ^              ********* 
                       ***.***          ********           *****.N*****
   ctr                 ***.***         ***********************...*****
    C - - - - - - - - -***E**D- - - - -Q****M...**********.....*****
                       ***.***         ********............******
                       ***.***           **********************
                      ***.***            ^   ***************
                     ***.***            ^ 
                    **A***             ^  
                     ****             ^    perp zone
                             zone 0  ^
                                    ^
                                   ^
*/
void np02_arc::get_distance_from_arc_zone_0(const np02_arc *a,
    np02_dist_from_xy_xy *result) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AA_ALWAYS_ASSERT(NULL != a);
AA_ALWAYS_ASSERT(NULL != result);
/* E = point on this arc centerline nearest to M = a->p_0 */
np02_xy xy_e(0.0,0.0);
const double dist_m_e = get_distance_from_xy_hw(a->m_p_0, 0.0, &xy_e);
//if( dist_m_e > 0.0 ) {
    const double dist_z_0 = a->distance_in_zone_0(xy_e);
    const double small_zone_error = get_small_distance() + (m_width/4);
    if (dist_z_0 >= -small_zone_error) {
        result->m_answer_quality = NP02_ANSWER_QUALITY_SMALL_ERROR;
        }
    else{
        if (dist_z_0 >= (-small_zone_error * 2.0)) {
            result->m_answer_quality = NP02_ANSWER_QUALITY_MEDIUM_ERROR;
            }
        else{
            result->m_answer_quality = NP02_ANSWER_QUALITY_LARGE_ERROR;            
            }        
        }
//    }
//else{
//    result->m_answer_quality = NP02_ANSWER_QUALITY_LARGE_ERROR;
//    }
const double hw = m_width / 2.0;
const double a_hw = (a->m_width)/2.0;
result->m_distance_from = dist_m_e - (hw + a_hw);
const np02_xy lambda_m_e = (a->m_p_0).get_unit_vector_to(xy_e);

/* D */
result->m_near_xy_defined = true;
(result->m_near_xy).set_x( xy_e.get_x() - (hw * lambda_m_e.get_x() ) );
(result->m_near_xy).set_y( xy_e.get_y() - (hw * lambda_m_e.get_y() ) );

/* Q */
result->m_other_near_xy_defined = true;
(result->m_other_near_xy).set_x( 
    (a->m_p_0).get_x() + (a_hw * lambda_m_e.get_x() ) );
(result->m_other_near_xy).set_y(
    (a->m_p_0).get_y() + (a_hw * lambda_m_e.get_y() ) );

AUTO_ASSERT( 0 == result->verify_result( this, a, 
    AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ) );
AA_DECR_CALL_DEPTH();
}

/*                                      
                                                   ^
                                                  ^
                     ****         zone 1         ^
                    **C***                      ^   perp
                     ***.***                   ^    zone       
                      ***.***              ** ^              
                       ***.***          ********           
   ctr                 ***.***         *************
    C - - - - - - - - -***E**D- - - - -Q****N.*********
                       ***.***         ********..*******
                       ***.***           ********...******
                      ***.***            ^   *******..******
                     ***.***            ^      ****** .******
                    **A***             ^         ******.******
                     ****             ^  perp     ******.******
                             zone 1  ^    zone     ******.******
                                    ^              ******.******
                                   P               ******.******
                                                   ******M******
                                                   *************
                                                    ***********
                                                        ***
*/
void np02_arc::get_distance_from_arc_zone_1(const np02_arc *a,
    np02_dist_from_xy_xy *result) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AA_ALWAYS_ASSERT(NULL != a);
AA_ALWAYS_ASSERT(NULL != result);
/* E = point on this arc centerline nearest to N = a->p_1 */
np02_xy xy_e(0.0,0.0);
const double dist_n_e = get_distance_from_xy_hw(a->m_p_1, 0.0, &xy_e);
//if( dist_n_e > 0.0 ) {
    const double dist_z_1 = a->distance_in_zone_1(xy_e);
    const double small_zone_error = get_small_distance() + (m_width/4);
    if (dist_z_1 >= -small_zone_error) {
        result->m_answer_quality = NP02_ANSWER_QUALITY_SMALL_ERROR;
        }
    else{
        const double small_distance = get_small_distance();
        if (dist_z_1 >= (-small_zone_error*2.0)) {
            result->m_answer_quality = NP02_ANSWER_QUALITY_MEDIUM_ERROR;
            }
        else{
            result->m_answer_quality = NP02_ANSWER_QUALITY_LARGE_ERROR;            
            }        
        }
//    }
//else{
//    result->m_answer_quality = NP02_ANSWER_QUALITY_LARGE_ERROR;
//    }
const double hw = m_width / 2.0;
const double a_hw = (a->m_width)/2.0;
result->m_distance_from = dist_n_e - (hw + a_hw);
const np02_xy lambda_n_e = (a->m_p_1).get_unit_vector_to(xy_e);

/* D */
result->m_near_xy_defined = true;
(result->m_near_xy).set_x( xy_e.get_x() - (hw * lambda_n_e.get_x() ) );
(result->m_near_xy).set_y( xy_e.get_y() - (hw * lambda_n_e.get_y() ) );

/* Q */
result->m_other_near_xy_defined = true;
(result->m_other_near_xy).set_x( 
    (a->m_p_1).get_x() + (a_hw * lambda_n_e.get_x() ) );
(result->m_other_near_xy).set_y(
    (a->m_p_1).get_y() + (a_hw * lambda_n_e.get_y() ) );

AUTO_ASSERT( 0 == result->verify_result( this, a, 
    AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ) );
AA_DECR_CALL_DEPTH();
}

/*                                      ^
                                          ^   
                                            ^  *****
                     ****          far      ********
                    **C***                 *****M****
                     ***.***             *****.*****
                      ***.***           *****.*****  ^
                       ***.***         *****.***** near^
   ctr                 ***.***         *****.*****       a->ctr
    C - - - - - - - - -***E**D- - - - -Q****.*****- - - - - P
                       ***.***         *****.***** near ^
                       ***.***  far    *****.*****  ^
                      ***.***           *****N*****
                     ***.***            ^*********
                    **A***          ^      *****
                     ****       ^
                            ^
inside perp_zone of a   ^
                    ^  outside perp_zone of a
*/
void np02_arc::get_distance_from_arc_far_perp_zone(const np02_arc *a,
    np02_dist_from_xy_xy *result) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AA_ALWAYS_ASSERT(NULL != a);
AA_ALWAYS_ASSERT(NULL != result);
/* E = point on this arc centerline nearest to P = a->ctr */
np02_xy xy_e(0.0,0.0);
const double dist_p_e = get_distance_from_xy_hw(a->m_ctr, 0.0, &xy_e);
if( ( a->is_in_perp_zone(xy_e) ) && ( dist_p_e > a->m_radius) && 
    ( !arc_is_approx_straight_line_seg() ) && 
    ( !(a->arc_is_approx_straight_line_seg()) ) ){
    result->m_answer_quality = NP02_ANSWER_QUALITY_SMALL_ERROR;
    }
else{
    result->m_answer_quality = NP02_ANSWER_QUALITY_LARGE_ERROR;
    }
const double hw = m_width / 2.0;
const double big_r_a = (a->m_radius) + ((a->m_width)/2.0);
result->m_distance_from = dist_p_e - (hw + big_r_a);
const np02_xy lambda_p_e = (a->m_ctr).get_unit_vector_to(xy_e);

/* D */
result->m_near_xy_defined = true;
(result->m_near_xy).set_x( xy_e.get_x() - (hw * lambda_p_e.get_x() ) );
(result->m_near_xy).set_y( xy_e.get_y() - (hw * lambda_p_e.get_y() ) );

/* Q */
result->m_other_near_xy_defined = true;
(result->m_other_near_xy).set_x( 
    (a->m_ctr).get_x() + (big_r_a * lambda_p_e.get_x() ) );
(result->m_other_near_xy).set_y(
    (a->m_ctr).get_y() + (big_r_a * lambda_p_e.get_y() ) );

AUTO_ASSERT( 0 == result->verify_result( this, a, 
    AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ) );
AA_DECR_CALL_DEPTH();
}

void np02_arc::get_distance_from_arc_centerline_intersect(const np02_arc *a,
    np02_dist_from_xy_xy *result_0, np02_dist_from_xy_xy *result_1) const{
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != a);
AA_ALWAYS_ASSERT(NULL != result_0);
AA_ALWAYS_ASSERT(NULL != result_1);

np02_xy xy_intsct_0; 
np02_xy xy_intsct_1;
const size_t ctrln_intsct_count = arc_arc_centerline_intersect( a,
    &xy_intsct_0, &xy_intsct_1);
size_t i = 0;
for( i = 0; i < 2; ++i ){
    const np02_xy& xy_intsct = ( 0 == i ) ? xy_intsct_0 : xy_intsct_1;
    np02_dist_from_xy_xy *result = ( 0 == i ) ? result_0 : result_1;
    if (i < ctrln_intsct_count){
        result->m_answer_quality = NP02_ANSWER_QUALITY_SMALL_ERROR;
        const double d_t_0_sq = xy_intsct.get_dsq_to(m_p_0);
        const double d_t_1_sq = xy_intsct.get_dsq_to(m_p_1);
        const double d_t=(d_t_0_sq<d_t_1_sq) ? sqrt(d_t_0_sq) : sqrt(d_t_1_sq);
        const double d_a_0_sq = xy_intsct.get_dsq_to(a->m_p_0);
        const double d_a_1_sq = xy_intsct.get_dsq_to(a->m_p_1);
        const double d_a=(d_a_0_sq<d_a_1_sq) ? sqrt(d_a_0_sq) : sqrt(d_a_1_sq);
        result->m_distance_from = -(d_t + d_a + ((m_width + a->m_width)/2.0));
        result->m_near_xy_defined = true;
        result->m_other_near_xy_defined = true;
        result->m_near_xy = xy_intsct;
        result->m_other_near_xy = xy_intsct;
        }
    else{
        result->m_answer_quality = NP02_ANSWER_QUALITY_LARGE_ERROR;
        result->m_distance_from = np02_shape::m_googol;
        result->m_near_xy_defined = false;
        result->m_other_near_xy_defined = false;
        (result->m_near_xy).set_xy(0.0,0.0);
        (result->m_other_near_xy).set_xy(0.0,0.0);
        }

    AUTO_ASSERT( 0 == result->verify_result( this, a, 
        AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ) );
    }

AA_DECR_CALL_DEPTH();
}

void np02_arc::get_distance_from_shape_arc_p01(const np02_shape *shp,
    const bool& is_from_p0, np02_dist_from_xy_xy *result) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AA_ALWAYS_ASSERT(NULL != shp);
AA_ALWAYS_ASSERT(NULL != result);
const np02_xy p01 = ( is_from_p0 ) ? m_p_0 : m_p_1;
const double d_shp_from_p01 =
    shp->get_distance_from_xy( p01, &(result->m_other_near_xy) );

#if 0
const double dist_z = ( is_from_p0 ) ? 
    distance_in_zone_0(result->m_other_near_xy) :
    distance_in_zone_1(result->m_other_near_xy);
if (dist_z >= -m_width/4.0 ) {
    result->m_answer_quality = NP02_ANSWER_QUALITY_SMALL_ERROR;
    }
else{
    const double small_distance = get_small_distance();
    if (dist_z >= -small_distance) {
        result->m_answer_quality = NP02_ANSWER_QUALITY_MEDIUM_ERROR;
        }
    else{
        result->m_answer_quality = NP02_ANSWER_QUALITY_LARGE_ERROR;            
        }        
    }
const double hw = m_width / 2.0;
result->m_distance_from = d_shp_from_p01 - hw;
result->m_near_xy_defined = true;
result->m_other_near_xy_defined = true;
const np02_xy lambda_pd = p01.get_unit_vector_to( result->m_other_near_xy );
(result->m_near_xy).set_x( p01.get_x() + ( hw * lambda_pd.get_x() ) );
(result->m_near_xy).set_y( p01.get_y() + ( hw * lambda_pd.get_y() ) );
#else
result->m_answer_quality = NP02_ANSWER_QUALITY_SMALL_ERROR;
const double hw = m_width / 2.0;
result->m_distance_from = d_shp_from_p01 - hw;
result->m_near_xy_defined = true;
result->m_other_near_xy_defined = true;
const np02_xy lambda_pd = p01.get_unit_vector_to( result->m_other_near_xy );
(result->m_near_xy).set_x( p01.get_x() + ( hw * lambda_pd.get_x() ) );
(result->m_near_xy).set_y( p01.get_y() + ( hw * lambda_pd.get_y() ) );
#endif

AUTO_ASSERT( 0 == result->verify_result( this, shp, 
    AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ) );
AA_DECR_CALL_DEPTH();
}

void np02_arc::get_distance_from_shape_far_perp_zone(const np02_shape *shp,
    np02_dist_from_xy_xy *result) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AA_ALWAYS_ASSERT(NULL != shp);
AA_ALWAYS_ASSERT(NULL != result);
const double d_shp_from_ctr =
    shp->get_distance_from_xy( m_ctr, &(result->m_other_near_xy) );
if( is_in_perp_zone(result->m_other_near_xy) && ( d_shp_from_ctr >= m_radius )
    && ( !arc_is_approx_straight_line_seg() ) ){
    result->m_answer_quality = NP02_ANSWER_QUALITY_SMALL_ERROR;
    }
else{
    result->m_answer_quality = NP02_ANSWER_QUALITY_LARGE_ERROR;
    }
const double hw = m_width / 2.0;
result->m_near_xy_defined = true;
result->m_other_near_xy_defined = true;
const np02_xy lambda_cd = m_ctr.get_unit_vector_to( result->m_other_near_xy );
if( d_shp_from_ctr >= m_radius ){
    const double big_r = m_radius + hw;
    result->m_distance_from = d_shp_from_ctr - big_r; 
    (result->m_near_xy).set_x( m_ctr.get_x() + ( big_r * lambda_cd.get_x() ) );
    (result->m_near_xy).set_y( m_ctr.get_y() + ( big_r * lambda_cd.get_y() ) );
    }
else{
    const double small_r = (m_radius > hw) ? ( m_radius - hw ) : 0.0;
    result->m_distance_from = small_r - d_shp_from_ctr;
    (result->m_near_xy).set_x( m_ctr.get_x() + ( small_r * lambda_cd.get_x() ) );
    (result->m_near_xy).set_y( m_ctr.get_y() + ( small_r * lambda_cd.get_y() ) );
    }

AUTO_ASSERT( 0 == result->verify_result( this, shp, 
    AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ) );
AA_DECR_CALL_DEPTH();
}

void np02_arc::get_distance_from_line_seg_far_perp_zone(const np02_line_seg *n,
    np02_dist_from_xy_xy *result) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AA_ALWAYS_ASSERT(NULL != n);
AA_ALWAYS_ASSERT(NULL != result);
np02_xy seg_ctr_near_xy(0.0,0.0);
const double d_seg_ctr_from_ctr =
    n->get_distance_from_xy_hw( m_ctr, 0.0, &seg_ctr_near_xy );
const double hw = m_width / 2.0;
const double tol = hw/256.0;
if( is_in_perp_zone(seg_ctr_near_xy) && ( (d_seg_ctr_from_ctr + tol) >= m_radius )
    && ( !arc_is_approx_straight_line_seg() ) ){
    result->m_answer_quality = NP02_ANSWER_QUALITY_SMALL_ERROR;
    }
else{
    result->m_answer_quality = NP02_ANSWER_QUALITY_LARGE_ERROR;
    }
const double big_r = m_radius + hw;
const double seg_hw = (n->get_width())/2.0;
const double d_seg_ctr_from_arc = d_seg_ctr_from_ctr - big_r;
result->m_distance_from = d_seg_ctr_from_arc - seg_hw;
result->m_near_xy_defined = true;
result->m_other_near_xy_defined = true;
const np02_xy lambda_cd = m_ctr.get_unit_vector_to( seg_ctr_near_xy );
if( d_seg_ctr_from_arc < 0.0 ){
    result->m_other_near_xy = seg_ctr_near_xy;
    }
else{
    (result->m_other_near_xy).set_x( seg_ctr_near_xy.get_x() - ( seg_hw * lambda_cd.get_x() ) ); 
    (result->m_other_near_xy).set_y( seg_ctr_near_xy.get_y() - ( seg_hw * lambda_cd.get_y() ) ); 
    }
(result->m_near_xy).set_x( m_ctr.get_x() + ( big_r * lambda_cd.get_x() ) );
(result->m_near_xy).set_y( m_ctr.get_y() + ( big_r * lambda_cd.get_y() ) );

AUTO_ASSERT( 0 == result->verify_result( this, n, 
    AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ) );
AA_DECR_CALL_DEPTH();
}

void np02_arc::get_distance_from_rect_pt(const np02_rect *rect,
    const np02_rect_pt_idx& rect_pt_idx, np02_dist_from_xy_xy *result) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AA_ALWAYS_ASSERT(NULL != rect);
AA_ALWAYS_ASSERT(NULL != result);
const np02_xy& p = rect->get_pt_by_idx( rect_pt_idx );
const double d_from_p = get_distance_from_xy( p, &(result->m_near_xy) );
if( NP02_RECT_PT_IDX_CTR == rect_pt_idx ){
    result->m_distance_from = rect->get_distance_from_xy(
        result->m_near_xy, &(result->m_other_near_xy));
    }
else{
    result->m_distance_from = d_from_p;
    result->m_other_near_xy = p;
    }
result->m_answer_quality = NP02_ANSWER_QUALITY_SMALL_ERROR;
result->m_near_xy_defined = true;
result->m_other_near_xy_defined = true;

AUTO_ASSERT( 0 == result->verify_result( this, rect, 
    AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ) );
AA_DECR_CALL_DEPTH();
}

void np02_arc::get_distance_from_line_seg_p01(const np02_line_seg *n,
    const bool& is_from_seg_p0, np02_dist_from_xy_xy *result) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AA_ALWAYS_ASSERT(NULL != n);
AA_ALWAYS_ASSERT(NULL != result);
np02_xy ctr_ln_near_xy;
const np02_xy& seg_p_01 = (is_from_seg_p0) ? n->get_p_0() : n->get_p_1();
const double d_ctr_ln_from_seg_p_01 =
    get_distance_from_xy_hw(seg_p_01, 0.0, &ctr_ln_near_xy);
const double seg_len = n->get_len();
AUTO_ASSERT(seg_len >= 0.0);
const double near_end_tol = seg_len / 4.0;
const double fwd_dot_near_xy = (n->get_fwd()).dot(ctr_ln_near_xy);
const double near_seg_zone_01 = (is_from_seg_p0) ? 
    ( fwd_dot_near_xy <= ( n->get_fwd_dot_0() + near_end_tol ) ) :
    ( ( fwd_dot_near_xy + near_end_tol ) >= n->get_fwd_dot_1() );
if( near_seg_zone_01 ){
    result->m_answer_quality = NP02_ANSWER_QUALITY_SMALL_ERROR;
    }
else{
    result->m_answer_quality = NP02_ANSWER_QUALITY_LARGE_ERROR;
    }
const double hw = m_width / 2.0;
const double seg_hw = (n->get_width()) / 2.0;
result->m_distance_from = d_ctr_ln_from_seg_p_01 - (hw + seg_hw);
result->m_near_xy_defined = true;
result->m_other_near_xy_defined = true;
const np02_xy lambda_pd = ctr_ln_near_xy.get_unit_vector_to( seg_p_01 );
(result->m_near_xy).set_x(ctr_ln_near_xy.get_x() + (hw*lambda_pd.get_x()));
(result->m_near_xy).set_y(ctr_ln_near_xy.get_y() + (hw*lambda_pd.get_y()));
(result->m_other_near_xy).set_x(seg_p_01.get_x() - (seg_hw*lambda_pd.get_x()));
(result->m_other_near_xy).set_y(seg_p_01.get_y() - (seg_hw*lambda_pd.get_y()));

AUTO_ASSERT( 0 == result->verify_result( this, n, 
    AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ) );
AA_DECR_CALL_DEPTH();
}

void np02_arc::get_distance_from_line_seg_centerline_intersect(
    const np02_line_seg *n, np02_dist_from_xy_xy *result_0,
    np02_dist_from_xy_xy *result_1) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AA_ALWAYS_ASSERT(NULL != n);
AA_ALWAYS_ASSERT(NULL != result_0);
AA_ALWAYS_ASSERT(NULL != result_1);
np02_xy xy_intsct_0; 
np02_xy xy_intsct_1;
const size_t ctrln_intsct_count = arc_seg_centerline_intersect( n,
    &xy_intsct_0, &xy_intsct_1);
size_t i = 0;
for( i = 0; i < 2; ++i ){
    const np02_xy& xy_intsct = ( 0 == i ) ? xy_intsct_0 : xy_intsct_1;
    np02_dist_from_xy_xy *result = ( 0 == i ) ? result_0 : result_1;
    if (i < ctrln_intsct_count){
        result->m_answer_quality = NP02_ANSWER_QUALITY_SMALL_ERROR;
        const double d_t_0_sq = xy_intsct.get_dsq_to(m_p_0);
        const double d_t_1_sq = xy_intsct.get_dsq_to(m_p_1);
        const double d_t=(d_t_0_sq<d_t_1_sq) ? sqrt(d_t_0_sq) : sqrt(d_t_1_sq);
        const double d_n_0_sq = xy_intsct.get_dsq_to(n->get_p_0());
        const double d_n_1_sq = xy_intsct.get_dsq_to(n->get_p_1());
        const double d_n=(d_n_0_sq<d_n_1_sq) ? sqrt(d_n_0_sq) : sqrt(d_n_1_sq);
        result->m_distance_from = -(d_t + d_n + ((m_width + n->get_width())/2.0));
        result->m_near_xy_defined = true;
        result->m_other_near_xy_defined = true;
        result->m_near_xy = xy_intsct;
        result->m_other_near_xy = xy_intsct;
        }
    else{
        result->m_answer_quality = NP02_ANSWER_QUALITY_LARGE_ERROR;
        result->m_distance_from = np02_shape::m_googol;
        result->m_near_xy_defined = false;
        result->m_other_near_xy_defined = false;
        (result->m_near_xy).set_xy(0.0,0.0);
        (result->m_other_near_xy).set_xy(0.0,0.0);
        }
    
    AUTO_ASSERT( 0 == result->verify_result( this, n, 
        AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ) );
    }
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

const np02_xy& np02_rect::get_pt_by_idx( const np02_rect_pt_idx& i ) const{
const np02_xy *pt = NULL;
switch (i) {
    case NP02_RECT_PT_IDX_P00:
        pt = &m_p00;
        break;
    case NP02_RECT_PT_IDX_P01:
        pt = &m_p01;
        break;
    case NP02_RECT_PT_IDX_P10:
        pt = &m_p10;
        break;
    case NP02_RECT_PT_IDX_P11:
        pt = &m_p11;
        break;
    case NP02_RECT_PT_IDX_CTR:
    case NP02_RECT_PT_IDX_COUNT:
    default:
        pt = &m_ctr;
        break;
    }
AA_ALWAYS_ASSERT(NULL != pt);
return *pt;
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
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AA_XDBG_ASSERT(0 == verify_data_num(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_2);
double d = 0.0;
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
CF01_HASH_CONSISTENCY_CHECK( this->hash() );

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
if( dab_len_sq > np02_shape::m_teensy_ratio){
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
    if( dsq_ab > np02_shape::m_teensy_ratio_sq){
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
                        np_d_double_check * np02_shape::m_small_ratio :
                        np02_shape::m_small_ratio;
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
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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
        AUTO_ASSERT(fabs(dd - fabs(d_ctr)) <=
            fabs((dd + fabs(d_ctr))*np02_shape::m_little_ratio) );
        np02_xy unit_v(1.0,0.0);
        if( dd < np02_shape::m_teensy_ratio){
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

double np02_rect::get_distance_from_arc(const np02_arc *a,
    np02_xy *near_xy, np02_xy *arc_near_xy) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AA_ALWAYS_ASSERT(NULL != a);
double d = a->get_distance_from_rect( this, arc_near_xy, near_xy );
AA_DECR_CALL_DEPTH();
return d;
}

double np02_rect::get_distance_from_line_seg(
    const np02_line_seg *n, np02_xy *near_xy,
    np02_xy *line_seg_near_xy)const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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
    AUTO_ASSERT(fabs(dd - fabs(d_ctr)) <=
        fabs((dd + fabs(d_ctr))*np02_shape::m_little_ratio) );
    np02_xy unit_v(1.0,0.0);
    if( dd < np02_shape::m_teensy_ratio){
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
AA_DECR_CALL_DEPTH();
return d;
}

double np02_rect::get_distance_from_rect(const np02_rect *r,
    np02_xy *near_xy, np02_xy *rect_near_xy) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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
AA_DECR_CALL_DEPTH();
return d;
}

double np02_rect::get_distance_from_polygon(const np02_polygon *p,
    np02_xy *near_xy, np02_xy *rect_near_xy) const{
/* TODO: implement polygon */
AA_ALWAYS_ASSERT(false);
return 0.0;
}

double np02_rect::get_distance_from_shape(const np02_shape *s,
    np02_xy *near_xy, np02_xy *other_near_xy) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
double d = 0.0;
if(NULL != s){
    const np02_circle *c = dynamic_cast<const np02_circle *>(s);
    const np02_arc *a = dynamic_cast<const np02_arc *>(s);
    const np02_line_seg *n = dynamic_cast<const np02_line_seg *>(s);
    const np02_rect *r = dynamic_cast<const np02_rect *>(s);
    const np02_polygon *p = dynamic_cast<const np02_polygon *>(s);
    const np02_spline *i = dynamic_cast<const np02_spline *>(s);
    if(NULL != c){ d = get_distance_from_circle(c, near_xy, other_near_xy); }
    else if(NULL!=a){ d=get_distance_from_arc(a, near_xy, other_near_xy);}
    else if(NULL!=n){ d=get_distance_from_line_seg(n, near_xy, other_near_xy);}
    else if(NULL!=r){ d=get_distance_from_rect(r, near_xy, other_near_xy); }
    else if(NULL!=p){ d=get_distance_from_polygon(p, near_xy, other_near_xy); }
    else if(NULL!=i){ AA_ALWAYS_ASSERT(false);/* spline not supported */}
    else{ AA_ALWAYS_ASSERT(false); }
    }
AA_DECR_CALL_DEPTH();
return d;
}

void np02_rect::translate_no_loc_grid(const np02_xy& dxy){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
const np02_xy ctr(m_ctr.get_x() + dxy.get_x(), m_ctr.get_y() + dxy.get_y());
init(ctr, m_w, m_h, m_rot_deg);
AA_DECR_CALL_DEPTH();
}

void np02_rect::rotate_no_loc_grid(const np02_xy& rot_ctr,
    const double& rot_deg){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
double cos_rot = 1.0;
double sin_rot = 0.0;
cos_sin_rot_deg( rot_deg, &cos_rot, &sin_rot );
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
AA_DECR_CALL_DEPTH();
}

uint64_t np02_rect::hash( const uint64_t& h_in ) const{
uint64_t h = np02_shape::hash( h_in );
h = m_ctr.hash( h );
#if defined( CF01_SUPPORT )
h = cf01_obj_hash( h, m_w );
h = cf01_obj_hash( h, m_h );
h = cf01_obj_hash( h, m_rot_deg );
#else
h += reinterpret_cast<cf01_uint64>(m_w);
h ^= ((h << 17) | (h >> 47));
h += reinterpret_cast<cf01_uint64>(m_h);
h ^= ((h << 27) | (h >> 37));
h += reinterpret_cast<cf01_uint64>(m_rot_deg);
h ^= ((h << 37) | (h >> 27));
#endif
h = m_fwd.hash( h );
#if defined( CF01_SUPPORT )
h = cf01_obj_hash( h, m_fwd_dot_ctr );
h = cf01_obj_hash( h, m_fwd_cross_ctr );
#else
h += reinterpret_cast<cf01_uint64>(m_fwd_dot_ctr);
h ^= ((h << 47) | (h >> 17));
h += reinterpret_cast<cf01_uint64>(m_fwd_cross_ctr);
h ^= ((h << 57) | (h >> 7));
#endif
h = m_p00.hash( h );
h = m_p01.hash( h );
h = m_p10.hash( h );
h = m_p11.hash( h );
return h;
}

int np02_rect::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const{
int err_cnt = 0;
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
if(fabs(m_fwd.get_x() - cos_rot) > np02_shape::m_small_ratio){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }
if(fabs(m_fwd.get_y() - sin_rot) > np02_shape::m_small_ratio){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
    }

const double fwd_len_sq = (m_fwd.get_x() * m_fwd.get_x()) +
                          (m_fwd.get_y() * m_fwd.get_y());
if(fabs(fwd_len_sq - 1.0) > np02_shape::m_small_ratio){
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
const double max_d_sq_err = (len_sq_criteria > 1.0) ?
    (len_sq_criteria * np02_shape::m_little_ratio_sq) :
    np02_shape::m_little_ratio_sq;
const double max_d_err = (len_sq_criteria > 1.0) ?
    sqrt(len_sq_criteria) * np02_shape::m_small_ratio :
    np02_shape::m_small_ratio;

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
        "w1_sq-w_sq=%g  max_d_sq_err=%g\n", __FILE__, __LINE__, __FUNCTION__,
        dx, dy, w1_sq, w_sq, w1_sq-w_sq, max_d_sq_err );
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
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
np02_shape::write_bmp_file(xy_min, pixel_num, color, bmp_file);
AA_DECR_CALL_DEPTH();
}

void np02_rect::write_dxf_file(const std::string& layer,
    const uint8_t& color, np02_dxf_file *dxf_file) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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
if(dsq < np02_shape::m_teensy_ratio_sq){
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
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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

/* hw=0 => get distance from line segment centerline */
double np02_line_seg::get_distance_from_xy(const np02_xy& xy,
    np02_xy *near_xy) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
const double hw = m_width / 2.0;
const double d = get_distance_from_xy_hw( xy, hw, near_xy );
AA_DECR_CALL_DEPTH();
return d;
}

double np02_line_seg::get_distance_from_xy_hw(const np02_xy& xy,
    const double& hw, np02_xy *near_xy) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AUTO_ASSERT(m_width >= 0.0);
AA_XDBG_ASSERT(0 == verify_data_num(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_2);

double dx, dy, dsq;
double d_ctr = 0.0;
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
AA_ALWAYS_ASSERT(false);
return 0.0;
}

double np02_line_seg::get_distance_from_circle(
    const np02_circle *c, np02_xy *near_xy,
    np02_xy *circle_near_xy) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AA_ALWAYS_ASSERT(NULL != c);
AUTO_ASSERT(c->get_radius() >= 0.0);
AUTO_ASSERT(m_width >= 0.0);
AA_XDBG_ASSERT(0 == verify_data_num(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), AA_DEBUG_LEVEL_2);
const double d = c->get_distance_from_line_seg(this, circle_near_xy, near_xy);
AA_DECR_CALL_DEPTH();
return d;
}

double np02_line_seg::get_distance_from_arc(const np02_arc *a,
    np02_xy *near_xy, np02_xy *arc_near_xy) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AA_ALWAYS_ASSERT(NULL != a);
double d = a->get_distance_from_line_seg( this, arc_near_xy, near_xy );
AA_DECR_CALL_DEPTH();
return d;
}

double np02_line_seg::get_distance_from_line_seg(
    const np02_line_seg *n, np02_xy *near_xy,
    np02_xy *line_seg_near_xy)const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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
const bool intersect_possible = 
                       ( ( ( cross_m0n1m1 <= 0.0) &&
                           ( cross_n1m1n0 <= 0.0) &&
                           ( cross_m1n0m0 <= 0.0) &&
                           ( cross_n0m0n1 <= 0.0) ) ||
                         ( ( cross_m0n1m1 >= 0.0) &&
                           ( cross_n1m1n0 >= 0.0) &&
                           ( cross_m1n0m0 >= 0.0) &&
                           ( cross_n0m0n1 >= 0.0) ) );
if(intersect_possible){
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
    np02_xy p(0.0,0.0);
    if( fabs(det) > np02_shape::m_tiny_ratio){
        p.set_x((( n_fwd_x * m_fwd_cross_01 ) - 
                     ( m_fwd_x * n_fwd_cross_01 ))/det );
        p.set_y((( n_fwd_y * m_fwd_cross_01 ) - 
                     ( m_fwd_y * n_fwd_cross_01 ))/det );
        }
    else{
        /* take the average of all endpoints */
        p.set_x((xm0 + xm1 + xn0 + xn1)/4.0 );
        p.set_y((ym0 + ym1 + yn0 + yn1)/4.0 );
        }

    const double& m_fwd_dot_p = m_fwd.dot(p);
    const double& n_fwd_dot_p = (n->get_fwd()).dot(p);

    if( ( m_fwd_dot_0 <= m_fwd_dot_p ) &&
        ( m_fwd_dot_p <= m_fwd_dot_1 ) &&
        ( n->m_fwd_dot_0 <= n_fwd_dot_p ) &&
        ( n_fwd_dot_p <= n->m_fwd_dot_1 ) ){
    
    
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
        if((len_m > np02_shape::m_teensy_ratio) &&
           (len_n > np02_shape::m_teensy_ratio)){
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
                nr_xy = p;
                oth_nr_xy = p;
                }
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
    const double max_allowed_err = ( 2 * d2_err_estimate ) +
        ( ( fabs(d) + fabs(d2) )*np02_shape::m_small_ratio );

    if( d > 0.0 ){

        if(fabs(d - d2) > max_allowed_err){
            static int printf_count = 0;
            if( printf_count < 50 ){
                ++printf_count;
                std::cout << "fabs(d - d2) > (2 * d2_err_estimate)\n";
                std::cout << "d = " << d << "\n";
                std::cout << "d2 = " << d2 << "\n";
                std::cout << "d2_err_estimate = " << d2_err_estimate << "\n";
                std::cout << "max_allowed_err = " << max_allowed_err << "\n";
                std::cout << "nr_xy = " << nr_xy.get_x() << "," << nr_xy.get_y() << "\n";
                std::cout << "near_xy2 = " << near_xy2.get_x() << "," << near_xy2.get_y() << "\n";
                std::cout << "oth_nr_xy = " << oth_nr_xy.get_x() << "," << oth_nr_xy.get_y() << "\n";
                std::cout << "line_seg_near_xy2 = " << line_seg_near_xy2.get_x()
                    << "," << line_seg_near_xy2.get_y() << "\n";
                }
            }

        AA_XDBG_ASSERT(fabs(d - d2) <= max_allowed_err,CF01_AA_DEBUG_LEVEL_1 );
        AA_XDBG_ASSERT(near_xy_err <= max_allowed_err,CF01_AA_DEBUG_LEVEL_1 );
        AA_XDBG_ASSERT(line_seg_near_xy_err <= max_allowed_err,CF01_AA_DEBUG_LEVEL_1 );
        }
    }

AA_DECR_CALL_DEPTH();
return d;
}


double np02_line_seg::get_distance_from_line_seg_double_check(
    const np02_line_seg *n, double *err_estimate,
    np02_xy *near_xy, np02_xy *line_seg_near_xy)const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
AA_ALWAYS_ASSERT(NULL != r);
const double d = r->get_distance_from_line_seg(this, rect_near_xy, near_xy);
AA_DECR_CALL_DEPTH();
return d;
}

double np02_line_seg::get_distance_from_polygon(
    const np02_polygon *p,  np02_xy *near_xy,
    np02_xy *rect_near_xy) const{
/* TODO: implement polygon */
AA_ALWAYS_ASSERT(false);
return 0.0;
}

double np02_line_seg::get_distance_from_shape(const np02_shape *s,
    np02_xy *near_xy, np02_xy *other_near_xy) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
double d = 0.0;
if(NULL != s){
    const np02_circle *c = dynamic_cast<const np02_circle *>(s);
    const np02_arc *a = dynamic_cast<const np02_arc *>(s);
    const np02_line_seg *n = dynamic_cast<const np02_line_seg *>(s);
    const np02_rect *r = dynamic_cast<const np02_rect *>(s);
    const np02_polygon *p = dynamic_cast<const np02_polygon *>(s);
    const np02_spline *i = dynamic_cast<const np02_spline *>(s);
    if(NULL != c){ d = get_distance_from_circle(c, near_xy, other_near_xy); }
    else if(NULL != a){ d = get_distance_from_arc(a, near_xy, other_near_xy); }
    else if(NULL != n){ d=get_distance_from_line_seg(n,near_xy,other_near_xy);}
    else if(NULL != r){ d = get_distance_from_rect(r, near_xy, other_near_xy);}
    else if(NULL != p){ d=get_distance_from_polygon(p,near_xy,other_near_xy);}
    else if(NULL != i){ AA_ALWAYS_ASSERT(false);/* spline not supported */}
    else{ AA_ALWAYS_ASSERT(false); }
    }
AA_DECR_CALL_DEPTH();
return d;
}

void np02_line_seg::translate_no_loc_grid(const np02_xy& dxy){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
const np02_xy p0(m_p_0.get_x() + dxy.get_x(), m_p_0.get_y() + dxy.get_y());
const np02_xy p1(m_p_1.get_x() + dxy.get_x(), m_p_1.get_y() + dxy.get_y());
AUTO_ASSERT( fabs( m_p_0.get_dsq_to( m_p_1 ) - p0.get_dsq_to( p1 ) )
    <= m_small_ratio * ( m_p_0.get_dsq_to( m_p_1 ) + p0.get_dsq_to( p1 ) ) );
init(p0, p1, m_width);
AA_DECR_CALL_DEPTH();
}

void np02_line_seg::rotate_no_loc_grid(const np02_xy& rot_ctr,
const double& rot_deg){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
double cos_rot = 1.0;
double sin_rot = 0.0;
cos_sin_rot_deg( rot_deg, &cos_rot, &sin_rot );
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
    rot_ctr.get_x() + rot_arm_final_0.get_x(),
    rot_ctr.get_y() + rot_arm_final_0.get_y());
const np02_xy p1(
    rot_ctr.get_x() + rot_arm_final_1.get_x(),
    rot_ctr.get_y() + rot_arm_final_1.get_y());
AUTO_ASSERT( fabs( m_p_0.get_dsq_to( m_p_1 ) - p0.get_dsq_to( p1 ) )
    <= m_small_ratio * ( m_p_0.get_dsq_to( m_p_1 ) + p0.get_dsq_to( p1 ) ) );
init(p0, p1, m_width);
AA_DECR_CALL_DEPTH();
}

size_t np02_line_seg::seg_seg_centerline_intersect(const np02_line_seg *n,
    np02_xy *xy_intsct_0, np02_xy *xy_intsct_1) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
assert( NULL != n );
size_t intsct_count = 0;
/*
Solve
[  -a.fwd.y     a.fwd.x  ][ intsct.x ]  =  [ a.fwd_cross_01 ] 
[  -b.fwd.y     b.fwd.x  ][ intsct.y ]     [ b.fwd_cross_01 ] 

det= a.fwd.x*b.fwd.y - a.fwd.y*b.fwd.x

[ intsct.x ]  =  (1/det)   [   b.fwd.x   -a.fwd.x  ][ a.fwd_cross_01 ] 
[ intsct.y ]               [   b.fwd.y   -a.fwd.y  ][ b.fwd_cross_01 ] 
*/
const double det = ( m_fwd.get_x() * (n->m_fwd).get_y() ) -
                   ( m_fwd.get_y() * (n->m_fwd).get_x() );
if( fabs(det) > np02_shape::m_small_ratio ){
    const np02_xy p( ( ((n->m_fwd).get_x() * m_fwd_cross_01) - 
                       (m_fwd.get_x() * (n->m_fwd_cross_01)) )/det,
                     ( ((n->m_fwd).get_y() * m_fwd_cross_01) - 
                       (m_fwd.get_y() * (n->m_fwd_cross_01)) )/det );
    const double fwd_dot_p = m_fwd.dot(p);
    const double n_fwd_dot_p = (n->m_fwd).dot(p);
    if( ( m_fwd_dot_0 <= fwd_dot_p ) &&
        ( fwd_dot_p <= m_fwd_dot_1 ) &&
        ( n->m_fwd_dot_0 <= n_fwd_dot_p ) &&
        ( n_fwd_dot_p <= n->m_fwd_dot_1 ) ){
        intsct_count = 1;
        if( NULL != xy_intsct_0 ){
            *xy_intsct_0 = p;
            }
        }
    }
else{
    /* parallel */
    typedef std::pair<double, np02_xy> dist_xy;
    dist_xy dist_xy_array[4];
    int g = 0;
    for(g = 0; g < 4; ++g) {
        dist_xy *dstxy = &(dist_xy_array[g]);
        switch( g ){
            default:
            case 0:
                dstxy->first = fabs(get_distance_from_xy_hw(n->m_p_0,0.0, NULL));
                dstxy->second = n->m_p_0;
                break;
            case 1:
                dstxy->first = fabs(get_distance_from_xy_hw(n->m_p_1,0.0, NULL));
                dstxy->second = n->m_p_1;
                break;
            case 2:
                dstxy->first = fabs(n->get_distance_from_xy_hw(m_p_0,0.0, NULL));
                dstxy->second = m_p_0;
                break;
            case 3:
                dstxy->first = fabs(n->get_distance_from_xy_hw(m_p_1,0.0, NULL));
                dstxy->second = m_p_1;
                break;   
            }
        }
    std::sort(&(dist_xy_array[0]), &(dist_xy_array[4]) );
    const dist_xy& dstxy0 = dist_xy_array[0];
    if( dstxy0.first < np02_shape::m_little_ratio ){
        ++intsct_count;
        if( NULL != xy_intsct_0 ){
            *xy_intsct_0 = dstxy0.second;
            }

        const dist_xy& dstxy1 = dist_xy_array[1];
        if( dstxy1.first < np02_shape::m_little_ratio ){
            ++intsct_count;
            if( NULL != xy_intsct_1 ){
                *xy_intsct_1 = dstxy1.second;
                }
            }
        }
    }

AUTO_ASSERT( 0 == verify_data_seg_seg_centerline_intersect_result( n,
    xy_intsct_0, xy_intsct_1, intsct_count,
    AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ) );
AA_DECR_CALL_DEPTH();
return intsct_count;
}

int np02_line_seg::verify_data_seg_seg_centerline_intersect_result(
    const np02_line_seg *n, np02_xy *xy_intsct_0, np02_xy *xy_intsct_1,
    const size_t& intsct_count, char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
int err_cnt = 0;


const double seg_len = (get_fwd_dot_1()) - (get_fwd_dot_0());
const double n_seg_len = (n->get_fwd_dot_1()) - (n->get_fwd_dot_0());
const double len_sum = seg_len + n_seg_len;
const double max_d_err = ( len_sum > 1.0 ) ? 
    ( len_sum * m_small_ratio ) : m_small_ratio;

if( ( intsct_count > 0 ) && ( NULL != xy_intsct_0 ) ){

    const double d0_from_seg_0 = get_distance_from_xy_hw(
        *xy_intsct_0, 0.0, NULL );
    if( fabs(d0_from_seg_0) > max_d_err ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d0from_seg_0=%g != 0.0\n", d0_from_seg_0 ); 
        }

    const double d0_from_seg_1 = n->get_distance_from_xy_hw(
        *xy_intsct_0, 0.0, NULL );
    if( fabs(d0_from_seg_1) > max_d_err ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d0_from_seg_1=%g != 0.0\n", d0_from_seg_1 ); 
        }
    }

if( ( intsct_count > 1 ) && ( NULL != xy_intsct_1 ) ){

    const double d0_from_seg_0 = get_distance_from_xy_hw(
        *xy_intsct_0, 0.0, NULL );
    if( fabs(d0_from_seg_0) > max_d_err ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d0from_seg_0=%g != 0.0\n", d0_from_seg_0 ); 
        }

    const double d0_from_seg_1 = n->get_distance_from_xy_hw(
        *xy_intsct_0, 0.0, NULL );
    if( fabs(d0_from_seg_1) > max_d_err ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "d0_from_seg_1=%g != 0.0\n", d0_from_seg_1 ); 
        }
    }

if( err_cnt > 0 ){
    const np02_xy& seg1_p0 = n->get_p_0();
    const np02_xy& seg1_p1 = n->get_p_1();
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "verify_data_arc_seg_centerline_intersect_result("
        " seg0_p0=[%g,%g], seg0_p1=[%g,%g],"
        " seg1_p0=[%g,%g], seg1_p1=[%g,%g],"
        " xy_intsct_0=[%g,%g], xy_intsct_1=[%g,%g], intsct_count=%i\n",
        m_p_0.get_x(), m_p_0.get_y(), m_p_1.get_x(), m_p_1.get_y(), 
        seg1_p0.get_x(), seg1_p0.get_y(), seg1_p1.get_x(), seg1_p1.get_y(),
        ( NULL == xy_intsct_0 ) ? 0.0 : xy_intsct_0->get_x(),
        ( NULL == xy_intsct_0 ) ? 0.0 : xy_intsct_0->get_y(),
        ( NULL == xy_intsct_1 ) ? 0.0 : xy_intsct_1->get_x(),
        ( NULL == xy_intsct_1 ) ? 0.0 : xy_intsct_1->get_y(),
        intsct_count );
    }

AA_DECR_CALL_DEPTH();
return err_cnt;
}

uint64_t np02_line_seg::hash( const uint64_t& h_in ) const{
uint64_t h = np02_shape::hash( h_in );
h = m_p_0.hash( h );
h = m_p_1.hash( h );
#if defined( CF01_SUPPORT )
h = cf01_obj_hash( h, m_width );
#else
h += reinterpret_cast<cf01_uint64>(m_width);
h ^= ((h << 17) | (h >> 47));
#endif
h = m_fwd.hash( h );
#if defined( CF01_SUPPORT )
h = cf01_obj_hash( h, m_fwd_dot_0 );
h = cf01_obj_hash( h, m_fwd_dot_1 );
h = cf01_obj_hash( h, m_fwd_cross_01 );
#else
h += reinterpret_cast<cf01_uint64>(m_fwd_dot_0);
h ^= ((h << 27) | (h >> 37));
h += reinterpret_cast<cf01_uint64>(m_fwd_dot_1);
h ^= ((h << 37) | (h >> 27));
h += reinterpret_cast<cf01_uint64>(m_fwd_cross_01);
h ^= ((h << 47) | (h >> 17));
#endif
return h;
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
const double max_d01_err = (d01 > 1.0) ?
    d01 * np02_shape::m_small_ratio : np02_shape::m_small_ratio;
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

if(d01 > get_small_distance()){
    const np02_xy fwd(dx01/d01, dy01/d01);
    const double fwd_dot_0 = (fwd.get_x() * m_p_0.get_x()) +
                             (fwd.get_y() * m_p_0.get_y());
    const double fwd_dot_1 = (fwd.get_x() * m_p_1.get_x()) +
                             (fwd.get_y() * m_p_1.get_y());
    const double fwd_cross_0 = (fwd.get_x() * m_p_0.get_y()) -
                               (fwd.get_y() * m_p_0.get_x());
    const double fwd_cross_1 = (fwd.get_x() * m_p_1.get_y()) -
                               (fwd.get_y() * m_p_1.get_x());
    double dot_cross_err_threshold = np02_shape::m_small_ratio * (
        fabs(fwd_dot_0) + fabs(fwd_dot_1) + 
        fabs(fwd_cross_0) + fabs(fwd_cross_1) +
        fabs(m_fwd_dot_0) + fabs(m_fwd_dot_1) + 
        fabs(m_fwd_cross_01) );
    if(dot_cross_err_threshold < np02_shape::m_small_ratio){
        dot_cross_err_threshold = np02_shape::m_small_ratio;
        }
    if(fabs(m_fwd.get_x()-fwd.get_x()) > np02_shape::m_small_ratio){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "data error: %s[%i] %s\n", __FILE__, __LINE__, __FUNCTION__);
        }
    if(fabs(m_fwd.get_y()-fwd.get_y()) > np02_shape::m_small_ratio){
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
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
np02_shape::write_bmp_file(xy_min, pixel_num, color, bmp_file);
AA_DECR_CALL_DEPTH();
}

void np02_line_seg::write_dxf_file(const std::string& layer,
    const uint8_t& color, np02_dxf_file *dxf_file) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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
np02_xy_vec_itr unique_end = std::unique(m_vertices.begin(), m_vertices.end());
m_vertices.erase(unique_end, m_vertices.end());
}

void np02_polygon::np02_polygon::get_bb(np02_xy *xy_min,
    np02_xy *xy_max) const{
/* TODO: implement polygon */
AA_ALWAYS_ASSERT(false);
}

void np02_polygon::get_loc_grid_indices_for_init(
        const np02_loc_grid_dim& loc_grid_dim,
        const double& extra_search_d, np02_uint16_pair_vec *index_vec)
        const{
/* TODO: implement polygon */
AA_ALWAYS_ASSERT(false);
}

double np02_polygon::get_distance_from_xy(const np02_xy& xy,
    np02_xy *near_xy) const{
/* TODO: implement polygon */
AA_ALWAYS_ASSERT(false);
return 0.0;
}

double np02_polygon::get_distance_from_line_seg_ab(const np02_xy& xy_a,
    const np02_xy& xy_b, np02_xy *near_xy,
    np02_xy *other_near_xy) const{
/* TODO: implement polygon */
AA_ALWAYS_ASSERT(false);
return 0.0;
}

double np02_polygon::get_distance_from_circle(
    const np02_circle *c,
    np02_xy *near_xy, np02_xy *circle_near_xy) const{
/* TODO: implement polygon */
AA_ALWAYS_ASSERT(false);
return 0.0;
}

double np02_polygon::get_distance_from_arc(const np02_arc *a,
    np02_xy *near_xy, np02_xy *arc_near_xy) const{
/* TODO: implement polygon */
AA_ALWAYS_ASSERT(false);
return 0.0;
}

double np02_polygon::get_distance_from_line_seg(
    const np02_line_seg *n, np02_xy *near_xy,
    np02_xy *line_seg_near_xy)const{
/* TODO: implement polygon */
AA_ALWAYS_ASSERT(false);
return 0.0;
}

double np02_polygon::get_distance_from_rect(const np02_rect *r,
    np02_xy *near_xy, np02_xy *rect_near_xy) const{
/* TODO: implement polygon */
AA_ALWAYS_ASSERT(false);
return 0.0;
}

double np02_polygon::get_distance_from_polygon(
    const np02_polygon *p, np02_xy *near_xy,
    np02_xy *rect_near_xy) const{
/* TODO: implement polygon */
AA_ALWAYS_ASSERT(false);
return 0.0;
}

double np02_polygon::get_distance_from_shape(const np02_shape *s,
    np02_xy *near_xy, np02_xy *other_near_xy) const{
/* TODO: implement polygon */
AA_ALWAYS_ASSERT(false);
return 0.0;
}

void np02_polygon::translate_no_loc_grid(const np02_xy& dxy){
/* TODO: implement polygon */
AA_ALWAYS_ASSERT(false);
}

void np02_polygon::rotate_no_loc_grid(const np02_xy& rot_ctr,  const double& rot_deg){
/* TODO: implement polygon */
AA_ALWAYS_ASSERT(false);
}

uint64_t np02_polygon::hash( const uint64_t& h_in ) const{
uint64_t h = np02_shape::hash( h_in );
np02_xy_vec_citr v_itr = m_vertices.begin();
for( ; v_itr != m_vertices.end(); ++v_itr ){
    h = v_itr->hash(h);
    }
return h;
}

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



/* TODO: implement polygon */
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
/* TODO: implement polygon */
AA_ALWAYS_ASSERT(false);
}




np02_spline::np02_spline(){
set_shape_type(NP02_SHAPE_TYPE_SPLINE);
}

np02_spline::~np02_spline(){}

void np02_spline::np02_spline::get_bb(np02_xy *xy_min,
    np02_xy *xy_max) const{
/* TODO: implement spline */
AA_ALWAYS_ASSERT(false);
}

void np02_spline::get_loc_grid_indices_for_init(
        const np02_loc_grid_dim& loc_grid_dim,
        const double& extra_search_d, np02_uint16_pair_vec *index_vec)
        const{
/* TODO: implement spline */
AA_ALWAYS_ASSERT(false);
}

double np02_spline::get_distance_from_xy(const np02_xy& xy,
    np02_xy *near_xy) const{
/* TODO: implement spline */
AA_ALWAYS_ASSERT(false);
return 0.0;
}

double np02_spline::get_distance_from_line_seg_ab(const np02_xy& xy_a,
    const np02_xy& xy_b, np02_xy *near_xy,
    np02_xy *other_near_xy) const{
/* TODO: implement spline */
AA_ALWAYS_ASSERT(false);
return 0.0;
}

double np02_spline::get_distance_from_circle(
    const np02_circle *c,
    np02_xy *near_xy, np02_xy *circle_near_xy) const{
/* TODO: implement spline */
AA_ALWAYS_ASSERT(false);
return 0.0;
}

double np02_spline::get_distance_from_arc(const np02_arc *a,
    np02_xy *near_xy, np02_xy *arc_near_xy) const{
/* TODO: implement spline */
AA_ALWAYS_ASSERT(false);
return 0.0;
}

double np02_spline::get_distance_from_line_seg(
    const np02_line_seg *n, np02_xy *near_xy,
    np02_xy *line_seg_near_xy)const{
/* TODO: implement spline */
AA_ALWAYS_ASSERT(false);
return 0.0;
}

double np02_spline::get_distance_from_rect(const np02_rect *r,
    np02_xy *near_xy, np02_xy *rect_near_xy) const{
/* TODO: implement spline */
AA_ALWAYS_ASSERT(false);
return 0.0;
}

double np02_spline::get_distance_from_polygon(
    const np02_polygon *p, np02_xy *near_xy,
    np02_xy *rect_near_xy) const{
/* TODO: implement spline */
AA_ALWAYS_ASSERT(false);
return 0.0;
}

double np02_spline::get_distance_from_shape(const np02_shape *s,
    np02_xy *near_xy, np02_xy *other_near_xy) const{
/* TODO: implement spline */
AA_ALWAYS_ASSERT(false);
return 0.0;
}

void np02_spline::translate_no_loc_grid(const np02_xy& dxy){
/* TODO: implement spline */
AA_ALWAYS_ASSERT(false);
}

void np02_spline::rotate_no_loc_grid(const np02_xy& rot_ctr,  const double& rot_deg){
/* TODO: implement spline */
AA_ALWAYS_ASSERT(false);
}

uint64_t np02_spline::hash( const uint64_t& h_in ) const{
uint64_t h = np02_shape::hash( h_in );
return h;
}

int np02_spline::verify_data( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;
err_cnt += np02_shape::verify_data(err_msg, err_msg_capacity, err_msg_pos);

if(NP02_SHAPE_TYPE_SPLINE != get_shape_type() ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "spline: this=%x  shape_type=%i\n", this, get_shape_type()); 
    }

const np02_shp_alloc *shp_alloc = get_shp_alloc();
//if((NULL != shp_alloc) && 
//    (this != shp_alloc->alloc_get_spline_by_idx(get_shape_idx()))){
//    ++err_cnt;
//    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
//        "spline: this=%x  != (shp_alloc=%x) spline_by_idx(idx=%i)=%x\n",
//        this, shp_alloc, get_shape_idx(),
//        shp_alloc->alloc_get_spline_by_idx(get_shape_idx())); 
//    }



/* TODO: implement spline */
return err_cnt;
}

std::ostream& np02_spline::ostream_output(std::ostream& os) const{
return os;
}

void np02_spline::write_bmp_file(const np02_xy& xy_min,
    const double& pixel_num, const np02_bmp_color& color,
    np02_bmp_file *bmp_file) const{
np02_shape::write_bmp_file(xy_min, pixel_num, color, bmp_file);
}

void np02_spline::write_dxf_file(const std::string& layer,
    const uint8_t& color, np02_dxf_file *dxf_file) const{
/* TODO: implement spline */
AA_ALWAYS_ASSERT(false);
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

uint64_t np02_loc_grid_dim::hash( const uint64_t& h_in ) const{
uint64_t h = h_in;
#if defined( CF01_SUPPORT )
h = cf01_obj_hash( h, m_w );
h = cf01_obj_hash( h, m_h );
h = cf01_obj_hash( h, m_x_min );
h = cf01_obj_hash( h, m_y_min );
h = cf01_obj_hash( h, m_sq_size );
#else
h += reinterpret_cast<cf01_uint64>(m_w);
h ^= ((h << 17) | (h >> 47));
h += reinterpret_cast<cf01_uint64>(m_h);
h ^= ((h << 27) | (h >> 37));
h += reinterpret_cast<cf01_uint64>(m_x_min);
h ^= ((h << 37) | (h >> 27));
h += reinterpret_cast<cf01_uint64>(m_y_min);
h ^= ((h << 47) | (h >> 17));
h += reinterpret_cast<cf01_uint64>(m_sq_size);
h ^= ((h << 57) | (h >> 7));
#endif
return h;
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
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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
AA_DECR_CALL_DEPTH();
}

np02_shp_alloc *np02_loc_grid_node::get_shp_alloc() const{
return (NULL == m_loc_grid) ? NULL : m_loc_grid->get_shp_alloc();
}

lyr_idx_type np02_loc_grid_node::get_lyr_idx() const{
return (NULL == m_loc_grid) ?
    NP02_LYR_INVALID_IDX : m_loc_grid->get_lyr_idx();
}

uint64_t np02_loc_grid_node::hash( const uint64_t& h_in ) const{
uint64_t h = h_in;
#if defined( CF01_SUPPORT )
h = cf01_obj_hash( h, m_alloc_idx );
if( NULL != m_owner ){
    h = cf01_obj_hash( h, static_cast<uint8_t>(1) );
    }
if( NULL != m_loc_grid ){
    h = cf01_obj_hash( h, static_cast<uint8_t>(1) );
    }
h = cf01_obj_hash( h, m_i );
h = cf01_obj_hash( h, m_j );
if( NULL != m_prev ){
    h = cf01_obj_hash( h, static_cast<uint8_t>(1) );
    }
if( NULL != m_next ){
    h = cf01_obj_hash( h, static_cast<uint8_t>(1) );
    }
if( NULL != m_s_prev ){
    h = cf01_obj_hash( h, static_cast<uint8_t>(1) );
    }
if( NULL != m_s_next ){
    h = cf01_obj_hash( h, static_cast<uint8_t>(1) );
    }
#else
h += reinterpret_cast<cf01_uint64>(m_alloc_idx);
h ^= ((h << 7) | (h >> 57));
if( NULL != m_owner ){
    h += 1;
    h ^= ((h << 17) | (h >> 47));
    }
if( NULL != m_loc_grid ){
    h += 1;
    h ^= ((h << 27) | (h >> 37));
    }
h += reinterpret_cast<cf01_uint64>(m_i);
h ^= ((h << 37) | (h >> 27));
h += reinterpret_cast<cf01_uint64>(m_j);
h ^= ((h << 47) | (h >> 17));
if( NULL != m_prev ){
    h += 1;
    h ^= ((h << 57) | (h >> 7));
    }
if( NULL != m_next ){
    h += 1;
    h ^= ((h << 7) | (h >> 57));
    }
if( NULL != m_s_prev ){
    h += 1;
    h ^= ((h << 17) | (h >> 47));
    }
if( NULL != m_s_next ){
    h += 1;
    h ^= ((h << 27) | (h >> 37));
    }
#endif
return h;
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
        max_loc_grid_node_count = loc_grid_shp_alloc->alloc_get_total_shape_count();
        }
    }

if(NULL != m_owner){
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
    m_idx_shape_vec(), m_idx_pair_vec(), m_small_distance(0.0)
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
size_t loc_grid_sz;
if(!m_loc_grid_vec.empty()){
    clear_loc_grid();}
AUTO_ASSERT(m_loc_grid_vec.empty());
m_loc_grid_dim = d;
loc_grid_sz = static_cast<size_t>(d.get_w()) * static_cast<size_t>(d.get_h());
m_loc_grid_vec.resize(loc_grid_sz, NULL);
const double& sq_sz = m_loc_grid_dim.get_sq_size();
m_small_distance = ( sq_sz > 0.0 ) ?
    ( np02_shape::m_small_ratio * sq_sz ) : np02_shape::m_small_ratio;
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}

void np02_loc_grid::insert_shape_in_loc_grid(np02_shape *shape){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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
    np02_shape_vec *shapes) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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
    const np02_xy& xy_max, np02_shape_vec *shapes ) const{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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

uint64_t np02_loc_grid::hash( const uint64_t& h_in ) const{
uint64_t h = h_in;
#if defined( CF01_SUPPORT )
if( NULL != m_shp_alloc ){
    h = cf01_obj_hash( h, static_cast<uint8_t>(1) );
    }
h = cf01_obj_hash( h, m_alloc_idx );
h = cf01_obj_hash( h, m_lyr_idx );
#else
if( NULL != m_shp_alloc ){
    h += 1;
    h ^= ((h << 17) | (h >> 47));
    }
h += reinterpret_cast<cf01_uint64>(m_alloc_idx);
h ^= ((h << 7) | (h >> 57));
h += reinterpret_cast<cf01_uint64>(m_lyr_idx);
h ^= ((h << 17) | (h >> 47));
#endif
h = m_loc_grid_dim.hash(h);
#if defined( CF01_SUPPORT )
h = cf01_obj_hash( h, m_extra_search_d );
#else
h += reinterpret_cast<cf01_uint64>(m_extra_search_d);
h ^= ((h << 7) | (h >> 57));
#endif
loc_grid_node_vec_citr lgn_itr = m_loc_grid_vec.begin();
for( ; lgn_itr != m_loc_grid_vec.end(); ++lgn_itr ){
    const np02_loc_grid_node *lg_node= *lgn_itr;
    const size_t lg_count_max = 1000000;
    size_t lg_count = 0;
    while((lg_count < lg_count_max) && (NULL != lg_node)){
        ++lg_count;
        h = lg_node->hash(h);
        lg_node = lg_node->get_next();
        }
    }
#if defined( CF01_SUPPORT )
h = cf01_obj_hash( h, m_small_distance );
#else
h += reinterpret_cast<cf01_uint64>(m_small_distance);
h ^= ((h << 7) | (h >> 57));
#endif
return h;
}

int np02_loc_grid::verify_data( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;

size_t count, loc_grid_node_vec_sz_check;
uint16_t w=0, h=0, i=0, j=0;
size_t gn_vec_idx=0;
const np02_shape *shape=NULL;
const np02_loc_grid_node *loc_grid_node=NULL;

w=m_loc_grid_dim.get_w();
h=m_loc_grid_dim.get_h();
err_cnt=m_loc_grid_dim.verify_data(err_msg,err_msg_capacity,err_msg_pos);
loc_grid_node_vec_sz_check = static_cast<size_t>(w) * static_cast<size_t>(h);
size_t max_expected_shape_count= 1000000;
const np02_shp_alloc *shp_alloc = get_shp_alloc();
if(NULL != shp_alloc){
    max_expected_shape_count = shp_alloc->alloc_get_total_shape_count();

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
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
write_bmp_file_grid_lines(xy_min, pixel_num, color, bmp_file);
write_bmp_file_grid_info(xy_min, pixel_num, color, bmp_file);
AA_DECR_CALL_DEPTH();
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
        sprintf(msg_buf, "i%i j%i s%llu", i,j,
            static_cast<unsigned long long>(shape_count));

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
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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
CF01_HASH_CONSISTENCY_CHECK( static_cast<cf01_uint64>( 0 ) ); /* infinite loop debug */
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

const size_t np02_shp_alloc::m_max_free_chain_sz = 0x7FFFFFFF;

np02_shp_alloc::np02_shp_alloc():
    /* allocator */
    m_alloc_shape_handle_vec(),
    m_shape_handle_free_chain(NULL),
    m_alloc_circle_vec(),
    m_circle_free_chain(NULL),
    m_alloc_arc_vec(),
    m_arc_free_chain(NULL),
    m_alloc_line_seg_vec(),
    m_line_seg_free_chain(NULL),
    m_alloc_rect_vec(),
    m_rect_free_chain(NULL),
    m_alloc_polygon_vec(),
    m_polygon_free_chain(NULL),
    m_alloc_spline_vec(),
    m_spline_free_chain(NULL),
    m_alloc_loc_grid_node_vec(),
    m_loc_grid_node_free_chain(NULL),
    m_alloc_loc_grid_vec(),
    m_alloc_boundary_seg_vec(),
    m_boundary_seg_free_chain(NULL),
    m_alloc_boundary_vec(),
    m_boundary_free_chain(NULL),
    m_alloc_region_vec(),
    m_region_free_chain(NULL)
{
AA_INCR_CALL_DEPTH();
AUTO_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()));
AA_DECR_CALL_DEPTH();
}

np02_shp_alloc::~np02_shp_alloc(){
alloc_delete_all();
}


np02_shape_handle *np02_shp_alloc::alloc_shape_handle(){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
np02_shape_handle *shape_handle = NULL;
if(NULL == m_shape_handle_free_chain){
    shape_handle = new np02_shape_handle();
    shape_handle->set_alloc_idx(static_cast<shape_idx_type>(m_alloc_shape_handle_vec.size()));
    shape_handle->set_shp_alloc(this);
    m_alloc_shape_handle_vec.push_back(shape_handle);
    }
else{
    shape_handle = m_shape_handle_free_chain;
    m_shape_handle_free_chain = m_shape_handle_free_chain->get_free_chain_next();
    shape_handle->set_free_chain_next(NULL);
    }
AUTO_ASSERT(0 == shape_handle->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_DECR_CALL_DEPTH();
return shape_handle;
}

void np02_shp_alloc::free_shape_handle(np02_shape_handle *shape_handle){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( static_cast<cf01_uint64>( 0 ) ); /* infinite loop debug */
AA_ALWAYS_ASSERT( NULL != shape_handle );
AA_ALWAYS_ASSERT( this == shape_handle->get_shp_alloc() );
shape_handle->destruct();
if(NULL != m_shape_handle_free_chain){
    shape_handle->set_free_chain_next(m_shape_handle_free_chain);
    }
m_shape_handle_free_chain = shape_handle;
shape_handle->set_owner_idx(NP02_SHP_HNDL_OWNER_INVALID_IDX);
AA_DECR_CALL_DEPTH();
}

void np02_shp_alloc::free_shape(np02_shape *shape)
{
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( static_cast<cf01_uint64>( 0 ) ); /* infinite loop debug */
if(NULL != shape){
    np02_circle *circle=NULL;
    np02_arc *arc = NULL;
    np02_line_seg *line_seg = NULL;
    np02_rect *rect = NULL;
    np02_polygon *polygon = NULL;
    np02_spline *spline = NULL;
    if( (circle = dynamic_cast<np02_circle *>(shape)) != NULL ){
        free_circle(circle);
        }
    else if((arc=dynamic_cast<np02_arc *>(shape)) != NULL){
        free_arc(arc);
        }
    else if((line_seg=dynamic_cast<np02_line_seg *>(shape)) != NULL){
        free_line_seg(line_seg);
        }
    else if((rect=dynamic_cast<np02_rect *>(shape)) != NULL){
        free_rect(rect);
        }
    else if((polygon=dynamic_cast<np02_polygon *>(shape)) != NULL){
        free_polygon(polygon);
        }
    else if((spline=dynamic_cast<np02_spline *>(spline)) != NULL){
        free_spline(spline);
        }
    else{ AA_ALWAYS_ASSERT(false); }
    }
AA_DECR_CALL_DEPTH();
}

np02_circle *np02_shp_alloc::alloc_circle(){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
np02_circle *circle = NULL;
if(NULL == m_circle_free_chain){
    circle = new np02_circle();
    circle->set_shape_idx(static_cast<shape_idx_type>(m_alloc_circle_vec.size()));
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
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( static_cast<cf01_uint64>( 0 ) ); /* infinite loop debug */
AA_ALWAYS_ASSERT( NULL != circle );
AA_ALWAYS_ASSERT( this == circle->get_shp_alloc() );
circle->destruct();
if(NULL != m_circle_free_chain){
    circle->set_free_chain_next(m_circle_free_chain);
    }
m_circle_free_chain = circle;

/* TODO: decide whether to keep the assert or keep the assignment*/
AA_ALWAYS_ASSERT( NULL == circle->get_shape_owner() );
circle->set_shape_owner(NULL);

AA_DECR_CALL_DEPTH();
}

np02_arc *np02_shp_alloc::alloc_arc(){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
np02_arc *arc = NULL;
if(NULL == m_arc_free_chain){
    arc = new np02_arc();
    arc->set_shape_idx(static_cast<shape_idx_type>(m_alloc_arc_vec.size()));
    arc->set_shp_alloc(this);
    m_alloc_arc_vec.push_back(arc);
    }
else{
    arc = m_arc_free_chain;
    m_arc_free_chain = m_arc_free_chain->get_free_chain_next();
    arc->set_free_chain_next(NULL);
    std::cout << arc << "\n";
    }
AUTO_ASSERT(0 == arc->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_DECR_CALL_DEPTH();
return arc;
}

void np02_shp_alloc::free_arc(np02_arc *arc){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( static_cast<cf01_uint64>( 0 ) ); /* infinite loop debug */
AA_ALWAYS_ASSERT( NULL != arc );
AA_ALWAYS_ASSERT( this == arc->get_shp_alloc() );
arc->destruct();
if(NULL != m_arc_free_chain){
    arc->set_free_chain_next(m_arc_free_chain);
    }
m_arc_free_chain = arc;

/* TODO: decide whether to keep the assert or keep the assignment*/
AA_ALWAYS_ASSERT( NULL == arc->get_shape_owner() );
arc->set_shape_owner(NULL);

AA_DECR_CALL_DEPTH();
}

np02_line_seg *np02_shp_alloc::alloc_line_seg(){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
np02_line_seg *line_seg = NULL;
if(NULL == m_line_seg_free_chain){
    line_seg = new np02_line_seg();
    line_seg->set_shape_idx(static_cast<shape_idx_type>(m_alloc_line_seg_vec.size()));
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
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( static_cast<cf01_uint64>( 0 ) ); /* infinite loop debug */
AA_ALWAYS_ASSERT( NULL != line_seg );
AA_ALWAYS_ASSERT( this == line_seg->get_shp_alloc() );
line_seg->destruct();
if(NULL != m_line_seg_free_chain){
    line_seg->set_free_chain_next(m_line_seg_free_chain);
    }
m_line_seg_free_chain = line_seg;

/* TODO: decide whether to keep the assert or keep the assignment*/
AA_ALWAYS_ASSERT( NULL == line_seg->get_shape_owner() );
line_seg->set_shape_owner(NULL);

AA_DECR_CALL_DEPTH();
}

np02_rect *np02_shp_alloc::alloc_rect(){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
np02_rect *rect = NULL;
if(NULL == m_rect_free_chain){
    rect = new np02_rect();
    rect->set_shape_idx(static_cast<shape_idx_type>(m_alloc_rect_vec.size()));
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
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( static_cast<cf01_uint64>( 0 ) ); /* infinite loop debug */
AA_ALWAYS_ASSERT( NULL != rect );
AA_ALWAYS_ASSERT( this == rect->get_shp_alloc() );
rect->destruct();
if(NULL != m_rect_free_chain){
    rect->set_free_chain_next(m_rect_free_chain);
    }
m_rect_free_chain = rect;

/* TODO: decide whether to keep the assert or keep the assignment*/
AA_ALWAYS_ASSERT( NULL == rect->get_shape_owner() );
rect->set_shape_owner(NULL);

AA_DECR_CALL_DEPTH();
}

np02_polygon *np02_shp_alloc::alloc_polygon(){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
np02_polygon *polygon = NULL;
if(NULL == m_polygon_free_chain){
    polygon = new np02_polygon();
    polygon->set_shape_idx(static_cast<shape_idx_type>(m_alloc_polygon_vec.size()));
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
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( static_cast<cf01_uint64>( 0 ) ); /* infinite loop debug */
AA_ALWAYS_ASSERT( NULL != polygon );
AA_ALWAYS_ASSERT( this == polygon->get_shp_alloc() );
polygon->destruct();
if(NULL != m_polygon_free_chain){
    polygon->set_free_chain_next(m_polygon_free_chain);
    }
m_polygon_free_chain = polygon;

/* TODO: decide whether to keep the assert or keep the assignment*/
AA_ALWAYS_ASSERT( NULL == polygon->get_shape_owner() );
polygon->set_shape_owner(NULL);

AA_DECR_CALL_DEPTH();
}

np02_spline *np02_shp_alloc::alloc_spline(){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
np02_spline *spline = NULL;
if(NULL == m_spline_free_chain){
    spline = new np02_spline();
    spline->set_shape_idx(static_cast<shape_idx_type>(m_alloc_spline_vec.size()));
    spline->set_shp_alloc(this);
    m_alloc_spline_vec.push_back(spline);
    }
else{
    spline = m_spline_free_chain;
    m_spline_free_chain = m_spline_free_chain->get_free_chain_next();
    spline->set_free_chain_next(NULL);
    }
AUTO_ASSERT(0 == spline->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_DECR_CALL_DEPTH();
return spline;
}

void np02_shp_alloc::free_spline(np02_spline *spline){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( static_cast<cf01_uint64>( 0 ) ); /* infinite loop debug */
AA_ALWAYS_ASSERT( NULL != spline );
AA_ALWAYS_ASSERT( this == spline->get_shp_alloc() );
spline->destruct();
if(NULL != m_spline_free_chain){
    spline->set_free_chain_next(m_spline_free_chain);
    }
m_spline_free_chain = spline;

/* TODO: decide whether to keep the assert or keep the assignment*/
AA_ALWAYS_ASSERT( NULL == spline->get_shape_owner() );
spline->set_shape_owner(NULL);
 
AA_DECR_CALL_DEPTH();
}

size_t np02_shp_alloc::alloc_get_total_shape_count() const{ 
const size_t total_shape_count = 
    alloc_get_circle_count() +
    alloc_get_arc_count() +
    alloc_get_line_seg_count() +
    alloc_get_rect_count() +
    alloc_get_polygon_count() +
    alloc_get_spline_count();
return total_shape_count;
}

np02_loc_grid_node *np02_shp_alloc::alloc_loc_grid_node(){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( static_cast<cf01_uint64>( 0 ) ); /* infinite loop debug */
/* TODO: call node->destruct() (might be named differently.) to free resources */
if(NULL != m_loc_grid_node_free_chain){
    loc_grid_node->set_free_chain_next(m_loc_grid_node_free_chain);
    }
m_loc_grid_node_free_chain = loc_grid_node;
loc_grid_node->set_owner(NULL);
AA_DECR_CALL_DEPTH();
}

np02_loc_grid *np02_shp_alloc::alloc_loc_grid(){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
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
AA_INCR_CALL_DEPTH();
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
AA_DECR_CALL_DEPTH();
}

np02_boundary_seg *np02_shp_alloc::alloc_boundary_seg(){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
np02_boundary_seg *boundary_seg = NULL;
if(NULL == m_boundary_seg_free_chain){
    boundary_seg = new np02_boundary_seg();
    boundary_seg->set_alloc_idx(static_cast<shape_idx_type>(m_alloc_boundary_seg_vec.size()));
    boundary_seg->set_shp_alloc(this);
    m_alloc_boundary_seg_vec.push_back(boundary_seg);
    }
else{
    boundary_seg = m_boundary_seg_free_chain;
    m_boundary_seg_free_chain = m_boundary_seg_free_chain->get_free_chain_next();
    boundary_seg->set_free_chain_next(NULL);
    std::cout << boundary_seg << "\n";
    }
AUTO_ASSERT(0 == boundary_seg->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_DECR_CALL_DEPTH();
return boundary_seg;
}

void np02_shp_alloc::free_boundary_seg(np02_boundary_seg *boundary_seg){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( static_cast<cf01_uint64>( 0 ) ); /* infinite loop debug */
AA_ALWAYS_ASSERT( NULL != boundary_seg );
AA_ALWAYS_ASSERT( this == boundary_seg->get_shp_alloc() );
boundary_seg->destruct();
if(NULL != m_boundary_seg_free_chain){
    boundary_seg->set_free_chain_next(m_boundary_seg_free_chain);
    }
m_boundary_seg_free_chain = boundary_seg;
boundary_seg->set_owner(NULL);
AA_DECR_CALL_DEPTH();
}

np02_boundary *np02_shp_alloc::alloc_boundary(){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
np02_boundary *boundary = NULL;
if(NULL == m_boundary_free_chain){
    boundary = new np02_boundary();
    boundary->set_alloc_idx(static_cast<shape_idx_type>(m_alloc_boundary_vec.size()));
    boundary->set_shp_alloc(this);
    m_alloc_boundary_vec.push_back(boundary);
    }
else{
    boundary = m_boundary_free_chain;
    m_boundary_free_chain = m_boundary_free_chain->get_free_chain_next();
    boundary->set_free_chain_next(NULL);
    std::cout << boundary << "\n";
    }
AUTO_ASSERT(0 == boundary->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_DECR_CALL_DEPTH();
return boundary;
}

void np02_shp_alloc::free_boundary(np02_boundary *boundary){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( static_cast<cf01_uint64>( 0 ) ); /* infinite loop debug */
AA_ALWAYS_ASSERT( NULL != boundary );
AA_ALWAYS_ASSERT( this == boundary->get_shp_alloc() );
boundary->destruct();
if(NULL != m_boundary_free_chain){
    boundary->set_free_chain_next(m_boundary_free_chain);
    }
m_boundary_free_chain = boundary;
boundary->set_owner(NULL);
AA_DECR_CALL_DEPTH();
}

np02_region *np02_shp_alloc::alloc_region(){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( this->hash() );
np02_region *region = NULL;
if(NULL == m_region_free_chain){
    region = new np02_region();
    region->set_alloc_idx(static_cast<shape_idx_type>(m_alloc_region_vec.size()));
    region->set_shp_alloc(this);
    m_alloc_region_vec.push_back(region);
    }
else{
    region = m_region_free_chain;
    m_region_free_chain = m_region_free_chain->get_free_chain_next();
    region->set_free_chain_next(NULL);
    std::cout << region << "\n";
    }
AUTO_ASSERT(0 == region->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_DECR_CALL_DEPTH();
return region;
}

void np02_shp_alloc::free_region(np02_region *region){
AA_INCR_CALL_DEPTH();
CF01_HASH_CONSISTENCY_CHECK( static_cast<cf01_uint64>( 0 ) ); /* infinite loop debug */
AA_ALWAYS_ASSERT( NULL != region );
AA_ALWAYS_ASSERT( this == region->get_shp_alloc() );
region->destruct();
if(NULL != m_region_free_chain){
    region->set_free_chain_next(m_region_free_chain);
    }
m_region_free_chain = region;
region->set_owner(NULL);
AA_DECR_CALL_DEPTH();
}

uint64_t np02_shp_alloc::hash( const uint64_t& h_in ) const{
uint64_t h = h_in;
size_t free_chain_sz = 0;

const np02_shape_handle *shape_handle = NULL;
shape_handle_vec::const_iterator shape_handle_itr = m_alloc_shape_handle_vec.begin();
for( ; shape_handle_itr != m_alloc_shape_handle_vec.end(); ++shape_handle_itr ){
    shape_handle = *shape_handle_itr;
    if( NULL == shape_handle ){
        #if defined( CF01_SUPPORT )
        h = cf01_obj_hash( h, static_cast<uint8_t>(1) );
        #else
        h += 1;
        h ^= ((h << 47) | (h >> 17));
        #endif
        }
    else{
        h = shape_handle->hash(h);
        }
    }

free_chain_sz = 0;
shape_handle = m_shape_handle_free_chain;
while( ( NULL != shape_handle ) && ( free_chain_sz < m_max_free_chain_sz ) ){
    h = shape_handle->hash(h);
    shape_handle = shape_handle->get_free_chain_next();
    ++free_chain_sz;
    }


const np02_circle *circle = NULL;
circle_vec::const_iterator circle_itr = m_alloc_circle_vec.begin();
for( ; circle_itr != m_alloc_circle_vec.end(); ++circle_itr ){
    circle = *circle_itr;
    if( NULL == circle ){
        #if defined( CF01_SUPPORT )
        h = cf01_obj_hash( h, static_cast<uint8_t>(1) );
        #else
        h += 1;
        h ^= ((h << 47) | (h >> 17));
        #endif
        }
    else{
        h = circle->hash(h);
        }
    }

free_chain_sz = 0;
circle = m_circle_free_chain;
while( ( NULL != circle ) && ( free_chain_sz < m_max_free_chain_sz ) ){
    h = circle->hash(h);
    circle = circle->get_free_chain_next();
    ++free_chain_sz;
    }


const np02_arc *arc = NULL;
arc_vec::const_iterator arc_itr = m_alloc_arc_vec.begin();
for( ; arc_itr != m_alloc_arc_vec.end(); ++arc_itr ){
    arc = *arc_itr;
    if( NULL == arc ){
        #if defined( CF01_SUPPORT )
        h = cf01_obj_hash( h, static_cast<uint8_t>(1) );
        #else
        h += 1;
        h ^= ((h << 47) | (h >> 17));
        #endif
        }
    else{
        h = arc->hash(h);
        }
    }

free_chain_sz = 0;
arc = m_arc_free_chain;
while( ( NULL != arc ) && ( free_chain_sz < m_max_free_chain_sz ) ){
    h = arc->hash(h);
    arc = arc->get_free_chain_next();
    ++free_chain_sz;
    }


const np02_line_seg *line_seg = NULL;
line_seg_vec::const_iterator line_seg_itr = m_alloc_line_seg_vec.begin();
for( ; line_seg_itr != m_alloc_line_seg_vec.end(); ++line_seg_itr ){
    line_seg = *line_seg_itr;
    if( NULL == line_seg ){
        #if defined( CF01_SUPPORT )
        h = cf01_obj_hash( h, static_cast<uint8_t>(1) );
        #else
        h += 1;
        h ^= ((h << 47) | (h >> 17));
        #endif
        }
    else{
        h = line_seg->hash(h);
        }
    }

free_chain_sz = 0;
line_seg = m_line_seg_free_chain;
while( ( NULL != line_seg ) && ( free_chain_sz < m_max_free_chain_sz ) ){
    h = line_seg->hash(h);
    line_seg = line_seg->get_free_chain_next();
    ++free_chain_sz;
    }


const np02_rect *rect = NULL;
rect_vec::const_iterator rect_itr = m_alloc_rect_vec.begin();
for( ; rect_itr != m_alloc_rect_vec.end(); ++rect_itr ){
    rect = *rect_itr;
    if( NULL == rect ){
        #if defined( CF01_SUPPORT )
        h = cf01_obj_hash( h, static_cast<uint8_t>(1) );
        #else
        h += 1;
        h ^= ((h << 47) | (h >> 17));
        #endif
        }
    else{
        h = rect->hash(h);
        }
    }

free_chain_sz = 0;
rect = m_rect_free_chain;
while( ( NULL != rect ) && ( free_chain_sz < m_max_free_chain_sz ) ){
    h = rect->hash(h);
    rect = rect->get_free_chain_next();
    ++free_chain_sz;
    }


const np02_polygon *polygon = NULL;
polygon_vec::const_iterator polygon_itr = m_alloc_polygon_vec.begin();
for( ; polygon_itr != m_alloc_polygon_vec.end(); ++polygon_itr ){
    polygon = *polygon_itr;
    if( NULL == polygon ){
        #if defined( CF01_SUPPORT )
        h = cf01_obj_hash( h, static_cast<uint8_t>(1) );
        #else
        h += 1;
        h ^= ((h << 47) | (h >> 17));
        #endif
        }
    else{
        h = polygon->hash(h);
        }
    }

free_chain_sz = 0;
polygon = m_polygon_free_chain;
while( ( NULL != polygon ) && ( free_chain_sz < m_max_free_chain_sz ) ){
    h = polygon->hash(h);
    polygon = polygon->get_free_chain_next();
    ++free_chain_sz;
    }


const np02_spline *spline = NULL;
spline_vec::const_iterator spline_itr = m_alloc_spline_vec.begin();
for( ; spline_itr != m_alloc_spline_vec.end(); ++spline_itr ){
    spline = *spline_itr;
    if( NULL == spline ){
        #if defined( CF01_SUPPORT )
        h = cf01_obj_hash( h, static_cast<uint8_t>(1) );
        #else
        h += 1;
        h ^= ((h << 47) | (h >> 17));
        #endif
        }
    else{
        h = spline->hash(h);
        }
    }

free_chain_sz = 0;
spline = m_spline_free_chain;
while( ( NULL != spline ) && ( free_chain_sz < m_max_free_chain_sz ) ){
    h = spline->hash(h);
    spline = spline->get_free_chain_next();
    ++free_chain_sz;
    }


const np02_loc_grid_node *lg_node = NULL;
loc_grid_node_vec::const_iterator lg_node_itr =
    m_alloc_loc_grid_node_vec.begin();
for( ; lg_node_itr != m_alloc_loc_grid_node_vec.end(); ++lg_node_itr ){
    lg_node = *lg_node_itr;
    if( NULL == lg_node ){
        #if defined( CF01_SUPPORT )
        h = cf01_obj_hash( h, static_cast<uint8_t>(1) );
        #else
        h += 1;
        h ^= ((h << 47) | (h >> 17));
        #endif
        }
    else{
        h = lg_node->hash(h);
        }
    }

free_chain_sz = 0;
lg_node = m_loc_grid_node_free_chain;
while( ( NULL != lg_node ) && ( free_chain_sz < m_max_free_chain_sz ) ){
    h = lg_node->hash(h);
    lg_node = lg_node->get_free_chain_next();
    ++free_chain_sz;
    }


const np02_loc_grid *loc_grid = NULL;
loc_grid_vec::const_iterator loc_grid_itr = m_alloc_loc_grid_vec.begin();
for( ; loc_grid_itr != m_alloc_loc_grid_vec.end(); ++loc_grid_itr ){
    loc_grid = *loc_grid_itr;
    if( NULL == loc_grid ){
        #if defined( CF01_SUPPORT )
        h = cf01_obj_hash( h, static_cast<uint8_t>(1) );
        #else
        h += 1;
        h ^= ((h << 47) | (h >> 17));
        #endif
        }
    else{
        h = loc_grid->hash(h);
        }
    }


const np02_boundary_seg *boundary_seg = NULL;
boundary_seg_vec::const_iterator boundary_seg_itr = m_alloc_boundary_seg_vec.begin();
for( ; boundary_seg_itr != m_alloc_boundary_seg_vec.end(); ++boundary_seg_itr ){
    boundary_seg = *boundary_seg_itr;
    if( NULL == boundary_seg ){
        #if defined( CF01_SUPPORT )
        h = cf01_obj_hash( h, static_cast<uint8_t>(1) );
        #else
        h += 1;
        h ^= ((h << 47) | (h >> 17));
        #endif
        }
    else{
        h = boundary_seg->hash(h);
        }
    }

free_chain_sz = 0;
boundary_seg = m_boundary_seg_free_chain;
while( ( NULL != boundary_seg ) && ( free_chain_sz < m_max_free_chain_sz ) ){
    h = boundary_seg->hash(h);
    boundary_seg = boundary_seg->get_free_chain_next();
    ++free_chain_sz;
    }


const np02_boundary *boundary = NULL;
boundary_vec::const_iterator boundary_itr = m_alloc_boundary_vec.begin();
for( ; boundary_itr != m_alloc_boundary_vec.end(); ++boundary_itr ){
    boundary = *boundary_itr;
    if( NULL == boundary ){
        #if defined( CF01_SUPPORT )
        h = cf01_obj_hash( h, static_cast<uint8_t>(1) );
        #else
        h += 1;
        h ^= ((h << 47) | (h >> 17));
        #endif
        }
    else{
        h = boundary->hash(h);
        }
    }

free_chain_sz = 0;
boundary = m_boundary_free_chain;
while( ( NULL != boundary ) && ( free_chain_sz < m_max_free_chain_sz ) ){
    h = boundary->hash(h);
    boundary = boundary->get_free_chain_next();
    ++free_chain_sz;
    }


const np02_region *region = NULL;
region_vec::const_iterator region_itr = m_alloc_region_vec.begin();
for( ; region_itr != m_alloc_region_vec.end(); ++region_itr ){
    region = *region_itr;
    if( NULL == region ){
        #if defined( CF01_SUPPORT )
        h = cf01_obj_hash( h, static_cast<uint8_t>(1) );
        #else
        h += 1;
        h ^= ((h << 47) | (h >> 17));
        #endif
        }
    else{
        h = region->hash(h);
        }
    }

free_chain_sz = 0;
region = m_region_free_chain;
while( ( NULL != region ) && ( free_chain_sz < m_max_free_chain_sz ) ){
    h = region->hash(h);
    region = region->get_free_chain_next();
    ++free_chain_sz;
    }


return h;
}

int np02_shp_alloc::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const{
int err_cnt = 0;
err_cnt += verify_data_alloc_shape_handle(err_msg, err_msg_capacity, err_msg_pos );
err_cnt += verify_data_alloc_circle(err_msg, err_msg_capacity, err_msg_pos );
err_cnt += verify_data_alloc_arc(err_msg, err_msg_capacity, err_msg_pos );
err_cnt += verify_data_alloc_line_seg(err_msg, err_msg_capacity, err_msg_pos);
err_cnt += verify_data_alloc_rect(err_msg, err_msg_capacity, err_msg_pos );
err_cnt += verify_data_alloc_polygon(err_msg, err_msg_capacity, err_msg_pos );
err_cnt += verify_data_alloc_spline(err_msg, err_msg_capacity, err_msg_pos );
err_cnt+=verify_data_alloc_loc_grid_node(err_msg,err_msg_capacity,err_msg_pos);
err_cnt += verify_data_alloc_loc_grid(err_msg, err_msg_capacity, err_msg_pos );
err_cnt += verify_data_alloc_boundary_seg(err_msg, err_msg_capacity, err_msg_pos );
err_cnt += verify_data_alloc_boundary(err_msg, err_msg_capacity, err_msg_pos );
err_cnt += verify_data_alloc_region(err_msg, err_msg_capacity, err_msg_pos );
return err_cnt;
}

std::ostream& np02_shp_alloc::ostream_output(std::ostream& os) const{
os << "<shp_alloc><this>" << std::hex << this << std::dec << "</this>\n";

/* allocator */
os << "</alloc_shape_handle_count=" << m_alloc_shape_handle_vec.size() << ">" 
    //<< "</shape_handle_free_chain_size=" << get_shape_handle_free_chain_size() << ">"
    << "\n";
os << "</alloc_circle_count=" << m_alloc_circle_vec.size() << ">" 
    //<< "</circle_free_chain_size=" << get_circle_free_chain_size() << ">"
    << "\n";
os << "</alloc_arc_count=" << m_alloc_arc_vec.size() << ">" 
    //<< "</arc_free_chain_size=" << get_arc_free_chain_size() << ">"
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
os << "</alloc_spline_count=" << m_alloc_spline_vec.size() << ">" 
    //<< "</spline_free_chain_size=" << get_spline_free_chain_size() << ">"
    << "\n";
os << "</alloc_loc_grid_node_count=" << m_alloc_loc_grid_node_vec.size() << ">" 
    //<< "</loc_grid_node_free_chain_size=" << get_loc_grid_node_free_chain_size() << ">"
    << "\n";
os << "</alloc_loc_grid_count=" << m_alloc_loc_grid_vec.size() << ">" << "\n";
os << "</alloc_boundary_seg_count=" << m_alloc_boundary_seg_vec.size() << ">" 
    //<< "</boundary_seg_free_chain_size=" << get_boundary_seg_free_chain_size() << ">"
    << "\n";
os << "</alloc_boundary_count=" << m_alloc_boundary_vec.size() << ">" 
    //<< "</boundary_free_chain_size=" << get_boundary_free_chain_size() << ">"
    << "\n";
os << "</alloc_region_count=" << m_alloc_region_vec.size() << ">" 
    //<< "</region_free_chain_size=" << get_region_free_chain_size() << ">"
    << "\n";
os << "</shp_alloc>\n";

return os;
}

/* destructor implementation.  Delete all objects. */
void np02_shp_alloc::alloc_delete_all(){
AA_INCR_CALL_DEPTH();
while (NULL != m_shape_handle_free_chain) {
    np02_shape_handle *shape_handle = m_shape_handle_free_chain;
    m_shape_handle_free_chain = m_shape_handle_free_chain->get_free_chain_next();
    shape_handle->set_free_chain_next(NULL);
    delete shape_handle;
    }
m_alloc_shape_handle_vec.clear();

while (NULL != m_circle_free_chain) {
    np02_circle *circle = m_circle_free_chain;
    m_circle_free_chain = m_circle_free_chain->get_free_chain_next();
    circle->set_free_chain_next(NULL);
    delete circle;
    }
m_alloc_circle_vec.clear();

while (NULL != m_arc_free_chain) {
    np02_arc *arc = m_arc_free_chain;
    m_arc_free_chain = m_arc_free_chain->get_free_chain_next();
    arc->set_free_chain_next(NULL);
    delete arc;
    }
m_alloc_arc_vec.clear();

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

while (NULL != m_spline_free_chain) {
    np02_spline *spline = m_spline_free_chain;
    m_spline_free_chain = m_spline_free_chain->get_free_chain_next();
    spline->set_free_chain_next(NULL);
    delete spline;
    }
m_alloc_spline_vec.clear();

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

while (NULL != m_boundary_seg_free_chain) {
    np02_boundary_seg *boundary_seg = m_boundary_seg_free_chain;
    m_boundary_seg_free_chain = m_boundary_seg_free_chain->get_free_chain_next();
    boundary_seg->set_free_chain_next(NULL);
    delete boundary_seg;
    }
m_alloc_boundary_seg_vec.clear();

while (NULL != m_boundary_free_chain) {
    np02_boundary *boundary = m_boundary_free_chain;
    m_boundary_free_chain = m_boundary_free_chain->get_free_chain_next();
    boundary->set_free_chain_next(NULL);
    delete boundary;
    }
m_alloc_boundary_vec.clear();

while (NULL != m_region_free_chain) {
    np02_region *region = m_region_free_chain;
    m_region_free_chain = m_region_free_chain->get_free_chain_next();
    region->set_free_chain_next(NULL);
    delete region;
    }
m_alloc_region_vec.clear();

AA_DECR_CALL_DEPTH();
}

int np02_shp_alloc::verify_data_alloc_shape_handle( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;
size_t i;
const np02_shape_handle *shape_handle = NULL;
for(i = 0; i < m_alloc_shape_handle_vec.size(); ++i){
    shape_handle = m_alloc_shape_handle_vec.at(i);
    if(NULL == shape_handle){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  m_alloc_shape_handle_vec.at(%i)=NULL\n",
            this, i ); 
        }
    else{
        if(shape_handle->get_alloc_idx() != i){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  shape_handle(%i)=%x->alloc_idx=%i\n",
                this, i, shape_handle, shape_handle->get_alloc_idx() );
            }
        const np02_shp_alloc *shape_handle_shp_alloc = shape_handle->get_shp_alloc();
        if((NULL != shape_handle_shp_alloc) && (this != shape_handle_shp_alloc)){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  != shape_handle(%i)=%x -> shp_alloc=%x\n",
                this, i, shape_handle, shape_handle_shp_alloc ); 
            }
        }
    }

i=0;
shape_handle=m_shape_handle_free_chain;
while((NULL != shape_handle) && (i < m_alloc_shape_handle_vec.size())){
    const size_t shape_handle_alloc_idx = shape_handle->get_alloc_idx();
    if(m_alloc_shape_handle_vec.size() <= shape_handle_alloc_idx){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain shape_handle %i "
            "alloc_idx=%i >= vec sz=%i\n", this, i, shape_handle_alloc_idx,
            m_alloc_shape_handle_vec.size() ); 
        }
    else if(m_alloc_shape_handle_vec.at(shape_handle_alloc_idx) != shape_handle){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain shape_handle %i =%x"
            " != vec.at(alloc_idx=%i) = %x\n", this, i, shape_handle, shape_handle_alloc_idx,
            m_alloc_shape_handle_vec.at(shape_handle_alloc_idx) ); 
        }
    
    ++i;
    shape_handle = shape_handle->get_free_chain_next();
    }
return err_cnt;
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

int np02_shp_alloc::verify_data_alloc_arc( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;
size_t i;
const np02_arc *arc = NULL;
for(i = 0; i < m_alloc_arc_vec.size(); ++i){
    arc = m_alloc_arc_vec.at(i);
    if(NULL == arc){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  m_alloc_arc_vec.at(%i)=NULL\n",
            this, i ); 
        }
    else{
        if(arc->get_shape_idx() != i){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  arc(%i)=%x->shape_idx=%i\n",
                this, i, arc, arc->get_shape_idx() );
            }
        const np02_shp_alloc *arc_shp_alloc = arc->get_shp_alloc();
        if((NULL != arc_shp_alloc) && (this != arc_shp_alloc)){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  != arc(%i)=%x -> shp_alloc=%x\n",
                this, i, arc, arc_shp_alloc ); 
            }
        }
    }

i=0;
arc=m_arc_free_chain;
while((NULL != arc) && (i < m_alloc_arc_vec.size())){
    const size_t arc_shp_idx = arc->get_shape_idx();
    if(m_alloc_arc_vec.size() <= arc_shp_idx){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain arc %i "
            "shp_idx=%i >= vec sz=%i\n", this, i, arc_shp_idx,
            m_alloc_arc_vec.size() ); 
        }
    else if(m_alloc_arc_vec.at(arc_shp_idx) != arc){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain arc %i =%x"
            " != vec.at(shp_idx=%i) = %x\n", this, i, arc, arc_shp_idx,
            m_alloc_arc_vec.at(arc_shp_idx) ); 
        }
    
    ++i;
    arc = arc->get_free_chain_next();
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

int np02_shp_alloc::verify_data_alloc_spline( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;
size_t i;
const np02_spline *spline = NULL;
for(i = 0; i < m_alloc_spline_vec.size(); ++i){
    spline = m_alloc_spline_vec.at(i);
    if(NULL == spline){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  m_alloc_spline_vec.at(%i)=NULL\n",
            this, i ); 
        }
    else{
        if(spline->get_shape_idx() != i){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  spline(%i)=%x->shape_idx=%i\n",
                this, i, spline, spline->get_shape_idx() );
            }
        const np02_shp_alloc *spline_shp_alloc = spline->get_shp_alloc();
        if((NULL != spline_shp_alloc) && (this != spline_shp_alloc)){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  != spline(%i)=%x -> shp_alloc=%x\n",
                this, i, spline, spline_shp_alloc ); 
            }
        }
    }

i=0;
spline=m_spline_free_chain;
while((NULL != spline) && (i < m_alloc_spline_vec.size())){
    const size_t spline_shp_idx = spline->get_shape_idx();
    if(m_alloc_spline_vec.size() <= spline_shp_idx){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain spline %i "
            "shp_idx=%i >= vec sz=%i\n", this, i, spline_shp_idx,
            m_alloc_spline_vec.size() ); 
        }
    else if(m_alloc_spline_vec.at(spline_shp_idx) != spline){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain spline %i =%x"
            " != vec.at(shp_idx=%i) = %x\n", this, i, spline, spline_shp_idx,
            m_alloc_spline_vec.at(spline_shp_idx) ); 
        }
    
    ++i;
    spline = spline->get_free_chain_next();
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

int np02_shp_alloc::verify_data_alloc_boundary_seg( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;
size_t i;
const np02_boundary_seg *boundary_seg = NULL;
for(i = 0; i < m_alloc_boundary_seg_vec.size(); ++i){
    boundary_seg = m_alloc_boundary_seg_vec.at(i);
    if(NULL == boundary_seg){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  m_alloc_boundary_seg_vec.at(%i)=NULL\n",
            this, i ); 
        }
    else{
        if(boundary_seg->get_alloc_idx() != i){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  boundary_seg(%i)=%x->alloc_idx=%i\n",
                this, i, boundary_seg, boundary_seg->get_alloc_idx() );
            }
        const np02_shp_alloc *boundary_seg_shp_alloc = boundary_seg->get_shp_alloc();
        if((NULL != boundary_seg_shp_alloc) && (this != boundary_seg_shp_alloc)){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  != boundary_seg(%i)=%x -> shp_alloc=%x\n",
                this, i, boundary_seg, boundary_seg_shp_alloc ); 
            }
        }
    }

i=0;
boundary_seg=m_boundary_seg_free_chain;
while((NULL != boundary_seg) && (i < m_alloc_boundary_seg_vec.size())){
    const size_t boundary_seg_alloc_idx = boundary_seg->get_alloc_idx();
    if(m_alloc_boundary_seg_vec.size() <= boundary_seg_alloc_idx){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain boundary_seg %i "
            "alloc_idx=%i >= vec sz=%i\n", this, i, boundary_seg_alloc_idx,
            m_alloc_boundary_seg_vec.size() ); 
        }
    else if(m_alloc_boundary_seg_vec.at(boundary_seg_alloc_idx) != boundary_seg){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain boundary_seg %i =%x"
            " != vec.at(alloc_idx=%i) = %x\n", this, i, boundary_seg, boundary_seg_alloc_idx,
            m_alloc_boundary_seg_vec.at(boundary_seg_alloc_idx) ); 
        }
    
    ++i;
    boundary_seg = boundary_seg->get_free_chain_next();
    }
return err_cnt;
}

int np02_shp_alloc::verify_data_alloc_boundary( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;
size_t i;
const np02_boundary *boundary = NULL;
for(i = 0; i < m_alloc_boundary_vec.size(); ++i){
    boundary = m_alloc_boundary_vec.at(i);
    if(NULL == boundary){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  m_alloc_boundary_vec.at(%i)=NULL\n",
            this, i ); 
        }
    else{
        if(boundary->get_alloc_idx() != i){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  boundary(%i)=%x->alloc_idx=%i\n",
                this, i, boundary, boundary->get_alloc_idx() );
            }
        const np02_shp_alloc *boundary_shp_alloc = boundary->get_shp_alloc();
        if((NULL != boundary_shp_alloc) && (this != boundary_shp_alloc)){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  != boundary(%i)=%x -> shp_alloc=%x\n",
                this, i, boundary, boundary_shp_alloc ); 
            }
        }
    }

i=0;
boundary=m_boundary_free_chain;
while((NULL != boundary) && (i < m_alloc_boundary_vec.size())){
    const size_t boundary_alloc_idx = boundary->get_alloc_idx();
    if(m_alloc_boundary_vec.size() <= boundary_alloc_idx){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain boundary %i "
            "alloc_idx=%i >= vec sz=%i\n", this, i, boundary_alloc_idx,
            m_alloc_boundary_vec.size() ); 
        }
    else if(m_alloc_boundary_vec.at(boundary_alloc_idx) != boundary){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain boundary %i =%x"
            " != vec.at(alloc_idx=%i) = %x\n", this, i, boundary, boundary_alloc_idx,
            m_alloc_boundary_vec.at(boundary_alloc_idx) ); 
        }
    
    ++i;
    boundary = boundary->get_free_chain_next();
    }
return err_cnt;
}

int np02_shp_alloc::verify_data_alloc_region( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;
size_t i;
const np02_region *region = NULL;
for(i = 0; i < m_alloc_region_vec.size(); ++i){
    region = m_alloc_region_vec.at(i);
    if(NULL == region){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  m_alloc_region_vec.at(%i)=NULL\n",
            this, i ); 
        }
    else{
        if(region->get_alloc_idx() != i){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  region(%i)=%x->alloc_idx=%i\n",
                this, i, region, region->get_alloc_idx() );
            }
        const np02_shp_alloc *region_shp_alloc = region->get_shp_alloc();
        if((NULL != region_shp_alloc) && (this != region_shp_alloc)){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "shp_alloc: this=%x  != region(%i)=%x -> shp_alloc=%x\n",
                this, i, region, region_shp_alloc ); 
            }
        }
    }

i=0;
region=m_region_free_chain;
while((NULL != region) && (i < m_alloc_region_vec.size())){
    const size_t region_alloc_idx = region->get_alloc_idx();
    if(m_alloc_region_vec.size() <= region_alloc_idx){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain region %i "
            "alloc_idx=%i >= vec sz=%i\n", this, i, region_alloc_idx,
            m_alloc_region_vec.size() ); 
        }
    else if(m_alloc_region_vec.at(region_alloc_idx) != region){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "shp_alloc: this=%x  free chain region %i =%x"
            " != vec.at(alloc_idx=%i) = %x\n", this, i, region, region_alloc_idx,
            m_alloc_region_vec.at(region_alloc_idx) ); 
        }
    
    ++i;
    region = region->get_free_chain_next();
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
    m_iteration(0), m_rand_uint32(4294967291ul), m_temp_err_msg_pos(0),
    m_test_err_msg_pos(0){
memset( m_temp_err_msg, 0, sizeof( m_temp_err_msg ) );
memset( m_test_err_msg, 0, sizeof( m_test_err_msg ) );
}

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
    np02_shape_vec *shapes, std::vector<np02_bmp_color> *colors ){
AA_INCR_CALL_DEPTH();
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
        CF01_HASH_CONSISTENCY_CHECK( cf01_obj_hash( np02_shape_vec_hash( *shapes ),
            ( NULL != shp_alloc ) ? shp_alloc->hash() : 0 ) );
        }
    }
AA_DECR_CALL_DEPTH();
}


np02_shape *np02_shape_test::make_rand_shape( np02_shp_alloc *shp_alloc,
    const np02_xy& shp_ctr, const double& basic_length,
    np02_bmp_color *color ){
np02_shape *shape = NULL;

double x0 = 0.0;
double y0 = 0.0;
double x1 = 0.0;
double y1 = 0.0;
double r = 0.0;
double w = 0.0;
double angle0 = 0.0;
double angle1 = 0.0;

advance_rand();
if( (get_rand_int()%64) < 32 ){
    /* grid points, integer-multiple dimensions */
    advance_rand();
    x0 = basic_length*static_cast<double>((get_rand_int()%25)-12)/16.0;
    advance_rand();
    y0 = basic_length*static_cast<double>((get_rand_int()%25)-12)/16.0;
    advance_rand();
    x1 = basic_length*static_cast<double>((get_rand_int()%25)-12)/16.0;
    advance_rand();
    y1 = basic_length*static_cast<double>((get_rand_int()%25)-12)/16.0;
    advance_rand();
    r = basic_length*static_cast<double>(get_rand_int()%17)/32.0;
    advance_rand();
    w = basic_length*static_cast<double>(get_rand_int()%17)/64.0;
    }
else{
    /* floating point */
    advance_rand();
    x0 = get_rand_dbl(-basic_length*0.75, basic_length*0.75);
    advance_rand();
    y0 = get_rand_dbl(-basic_length*0.75, basic_length*0.75);
    advance_rand();
    x1 = get_rand_dbl(-basic_length*0.75, basic_length*0.75);
    advance_rand();
    y1 = get_rand_dbl(-basic_length*0.75, basic_length*0.75);
    advance_rand();
    r = get_rand_dbl(0.0, basic_length/2.0);
    advance_rand();
    w = get_rand_dbl(0.0, basic_length/4.0);
    }

advance_rand();
if( (get_rand_int()%64) < 32 ){
    /* grid points, integer-multiple dimensions */
    advance_rand();
    angle0 = static_cast<double>(get_rand_int()%24)*15.0;
    advance_rand();
    angle1 = static_cast<double>(get_rand_int()%24)*15.0;
    }
else{
    /* floating point */
    advance_rand();
    angle0 = get_rand_dbl(0.0, 360.0);
    advance_rand();
    angle1 = get_rand_dbl(0.0, 360.0);
    }


np02_bmp_color clr(0,0,0);
switch(get_rand_int()%13){
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
            fabs(x1), fabs(y1), angle0);
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
            np02_xy(shp_ctr.get_x()+x1,shp_ctr.get_y()+y1), w);
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
    case 9:
    case 10:
    case 11:
    case 12:
        {
        np02_arc *arc = (NULL == shp_alloc) ? 
            new np02_arc() : shp_alloc->alloc_arc();
        advance_rand();
        if( (get_rand_int()%64) < 32 ){
            np02_arc::init_params arc_init_params;
            arc_init_params.m_ctr = shp_ctr;
            arc_init_params.m_radius = r; 
            arc_init_params.m_start_angle_deg = angle0; 
            arc_init_params.m_end_angle_deg = angle1;
            arc_init_params.m_width = ( w/2.0 > r ) ? (r / 2.0) : w;
            arc->init(arc_init_params);
            }
        else{
            np02_arc::init3pt_params arc_init3pt_params;
            arc_init3pt_params.m_pa.set_x(shp_ctr.get_x()+x0);
            arc_init3pt_params.m_pa.set_y(shp_ctr.get_y()+y0);
            arc_init3pt_params.m_pb = shp_ctr;
            arc_init3pt_params.m_pc.set_x(shp_ctr.get_x()+x1);
            arc_init3pt_params.m_pc.set_y(shp_ctr.get_y()+y1);
            arc_init3pt_params.m_width = w;
            arc_init3pt_params.m_max_radius = 1000000.0 * basic_length;
            np02_arc::init_params arc_init_params;
            np02_arc::init3pt_aux_params arc_init3pt_aux_params;
            np02_arc::init3pt_to_init_params( arc_init3pt_params,
                &arc_init_params, &arc_init3pt_aux_params );
            arc->init_force_p01( arc_init_params, arc_init3pt_aux_params );
            }
        shape = arc;

        /* create random color */
        advance_rand();

        switch(get_rand_int()%3){
            default:
            case 0: clr=np02_bmp_color::bl_navy();break;
            case 1: clr=np02_bmp_color::bl_fuchsia();break;
            case 2: clr=np02_bmp_color::bl_purple ();break;
            }
        }
        break;
    }

if(NULL != color){ *color = clr; }
assert( NULL != shape );
return shape;
}

void np02_shape_test::free_shapes( np02_shp_alloc *shp_alloc,
    np02_shape_vec *shapes ){
int e = 0;
if(NULL != shp_alloc){
    e = shp_alloc->verify_data(AA_ERR_BUF(),
        AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR());
    AA_ALWAYS_ASSERT(0 == e);
    }

if (NULL != shapes){
    np02_shape_vec_citr shp_itr = shapes->begin();
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

void np02_shape_test::draw_debug_shapes( const bmp_debug_file_params& p){

if( NULL != p.m_bmp_file ){
    size_t i;
    size_t j;
    std::ostringstream debug_ss;

    std::vector<np02_bmp_color> color_vec;
    if(NULL != p.m_color_vec){
        for( i = 0; i < p.m_color_vec->size(); ++i){
            np02_bmp_color color = p.m_color_vec->at(i);
            if(p.m_err_xy_found){
                color.m_red /= 4;
                color.m_green /= 4;
                color.m_blue /= 4;
                }
            color_vec.push_back(color);
            }        
        }

    if( ( NULL != p.m_shapes ) && ( color_vec.size() >= p.m_shapes->size() ) ){
        for( i = 0; i < p.m_shapes->size(); ++i){
            const np02_shape *shape = p.m_shapes->at(i);
            if( ( NULL != shape ) && ( p.m_err_shape_a != shape ) &&
                ( p.m_err_shape_b != shape ) ){
                np02_bmp_color color = color_vec.at(i);
                color.m_blend = true;
                shape->write_bmp_file(p.m_xy_min, p.m_pixel_num, color, p.m_bmp_file);
                }
            }
        }
    if( NULL != p.m_err_shape_a ){
        np02_bmp_color color_a = np02_bmp_color::red();
        color_a.m_blend = true;
        p.m_err_shape_a->write_bmp_file(p.m_xy_min, p.m_pixel_num, color_a, p.m_bmp_file);
        }
    if( NULL != p.m_err_shape_b ){
        np02_bmp_color color_b = np02_bmp_color::blue();
        color_b.m_blend = true;
        p.m_err_shape_b->write_bmp_file(p.m_xy_min, p.m_pixel_num, color_b, p.m_bmp_file);
        }

    if( ( NULL != p.m_shapes ) &&
        ( p.m_draw_overlap_circle || p.m_draw_connecting_line ) ){
        for( i = 0; i < p.m_shapes->size(); ++i){
            const np02_shape *shape_i = p.m_shapes->at(i);
            for( j = i+1; j < p.m_shapes->size(); ++j){
                const np02_shape *shape_j = p.m_shapes->at(j);
                np02_xy xy_a, xy_b;
        
                const double d = shape_i->get_distance_from_shape(shape_j,
                    &xy_a, &xy_b );
        
                const int32_t a_ii = static_cast<int32_t>(
                    (xy_a.get_x() - p.m_xy_min.get_x()) * p.m_pixel_num);
                const int32_t a_jj = static_cast<int32_t>(
                    (xy_a.get_y() - p.m_xy_min.get_y()) * p.m_pixel_num);
                const int32_t b_ii = static_cast<int32_t>(
                    (xy_b.get_x() - p.m_xy_min.get_x()) * p.m_pixel_num);
                const int32_t b_jj = static_cast<int32_t>(
                    (xy_b.get_y() - p.m_xy_min.get_y()) * p.m_pixel_num);
    
                /* dark red circle for overlap*/
                if( p.m_draw_overlap_circle && ( d < 0.0 ) ){
                    np02_xy xy_ab_ave( (xy_a.get_x() + xy_b.get_x())/2.0,
                        (xy_a.get_y() + xy_b.get_y())/2.0 );
                    const double ii_ctr =
                        ( xy_ab_ave.get_x() - p.m_xy_min.get_x() ) * p.m_pixel_num;
                    const double jj_ctr =
                        ( xy_ab_ave.get_y() - p.m_xy_min.get_y() ) * p.m_pixel_num;
                    const double rr = fabs(d) * p.m_pixel_num;
        
                    p.m_bmp_file->draw_circle(ii_ctr,jj_ctr,rr,np02_bmp_color(64,0,0, true));
                    }
    
                /* line connecting nearby shapes*/
                if( p.m_draw_connecting_line &&
                    (fabs(d) < (1.25 * p.m_basic_length)) ){
                    p.m_bmp_file->draw_line(a_ii, a_jj, b_ii, b_jj,
                        np02_bmp_color(255,255,255,true));
                    }
                }
            }
        }

    if(p.m_err_xy_found){
        debug_ss << __FILE__ << "[" << __LINE__ << "] "
            << "shape error  itr=" << m_iteration
            << "  drawing=" << p.m_footer << "/" << p.m_header;

        if( ( NULL != p.m_err_xy_a ) && ( NULL != p.m_err_xy_b ) ){
            const np02_xy text_pos((p.m_err_xy_a->get_x() + p.m_err_xy_b->get_x())/2.0,
                                       (p.m_err_xy_a->get_y() + p.m_err_xy_b->get_y())/2.0);
            const double d_err_ab = p.m_err_xy_a->get_distance_to(*p.m_err_xy_b);
            double hh = d_err_ab * p.m_pixel_num;
            if(hh < 10){ hh = 10.0; }
            const double ii_text =
                ( text_pos.get_x() - p.m_xy_min.get_x() ) * p.m_pixel_num;
            const double jj_text =
                ( text_pos.get_y() - p.m_xy_min.get_y() ) * p.m_pixel_num;
            debug_ss << "  xy(" << text_pos.get_x() << "," << text_pos.get_y() << ")"
                << "  A(" << p.m_err_xy_a->get_x() << "," << p.m_err_xy_a->get_y() << ")"
                << "  B(" << p.m_err_xy_b->get_x() << "," << p.m_err_xy_b->get_y() << ")"
                << "  pixel(" << ii_text << "," << jj_text << ")\n";

            p.m_bmp_file->draw_text("Error", ii_text, jj_text, hh, 15.0, np02_bmp_color::white());
            }

        debug_ss << "SHAPE_A:\n";
        if( NULL == p.m_err_shape_a ) {
            debug_ss << "NULL\n";
            }
        else{
            p.m_err_shape_a->ostream_output( debug_ss );
            }
        debug_ss << "SHAPE_B:\n";
        if( NULL == p.m_err_shape_b ) {
            debug_ss << "NULL\n";
            }
        else{
            p.m_err_shape_b->ostream_output( debug_ss );
            }

        if( NULL != p.m_err_xy_a ){
            const int32_t err_ii_a = static_cast<int32_t>(floor(
                ( p.m_err_xy_a->get_x() - p.m_xy_min.get_x() ) * p.m_pixel_num));
            const int32_t err_jj_a = static_cast<int32_t>(floor(
                ( p.m_err_xy_a->get_y() - p.m_xy_min.get_y() ) * p.m_pixel_num));
            p.m_bmp_file->draw_cross( err_ii_a, err_jj_a, 9, np02_bmp_color::yellow() );
            }

        if( NULL != p.m_err_xy_b ){
            const int32_t err_ii_b = static_cast<int32_t>(floor(
                ( p.m_err_xy_b->get_x() - p.m_xy_min.get_x() ) * p.m_pixel_num));
            const int32_t err_jj_b = static_cast<int32_t>(floor(
                ( p.m_err_xy_b->get_y() - p.m_xy_min.get_y() ) * p.m_pixel_num));
            p.m_bmp_file->draw_x( err_ii_b, err_jj_b, 7, np02_bmp_color::silver() );
            }
        }

    p.m_bmp_file->draw_text(p.m_footer.c_str(), 25.0, (p.m_bmp_file->get_height_px())-75.0,
        50.0, 0.0, np02_bmp_color::gray());
    p.m_bmp_file->draw_text(p.m_header.c_str(), 25.0, 25.0, 50.0, 0.0,
        np02_bmp_color::gray());

    if( NULL != p.m_debug_str_out){
        p.m_debug_str_out->append( debug_ss.str() );
        }
    }
}

void np02_shape_test::print_temp_buf_to_test_buf(
    const char *filename, const int line_num ){
/* copy to test buffer */
np02_snprintf(&(m_test_err_msg[0]), TEST_ERR_MSG_CAPACITY, &m_test_err_msg_pos,
    "\n[itr=%i]  %s[%i]\n%s\n", m_iteration, filename, line_num,
    &(m_temp_err_msg[0])); 

/* clear temp buffer */
memset( m_temp_err_msg, 0, sizeof( m_temp_err_msg ) );
m_temp_err_msg_pos = 0;
}

void np02_shape_test::print_and_clear_test_buf(){
/* print test buffer */
std::cout << &(m_test_err_msg[0]);

/* clear test buffer */
memset( m_test_err_msg, 0, sizeof( m_test_err_msg ) );
m_test_err_msg_pos = 0;
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
if( err_cnt > 0 ){ print_and_clear_test_buf(); }
if( NULL != AA_ERR_BUF_POS_PTR() && *(AA_ERR_BUF_POS_PTR()) > 0 ){ 
    std::cout << "\n\nAA_ERR_BUF:" << AA_ERR_BUF();
    }
AA_DECR_CALL_DEPTH();
return err_cnt;
}

int np02_shape_test::execute_test_1_0(){
AA_INCR_CALL_DEPTH();
int err_cnt = 0;
int e = 0;

np02_shape_vec shape_vec;
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


np02_arc arc0;
np02_arc::init_params arc0_init_params;
arc0_init_params.m_ctr.set_x(5.0);
arc0_init_params.m_ctr.set_y(-5.0);
arc0_init_params.m_radius = 6.0; 
arc0_init_params.m_start_angle_deg=30.0; 
arc0_init_params.m_end_angle_deg=135.0;
arc0_init_params.m_width=2.0;
arc0.init(arc0_init_params);
shape_vec.push_back(&arc0);
color_vec.push_back(np02_bmp_color(0, 255, 0, false));
e = arc0.verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR());
AA_ALWAYS_ASSERT(0 == e);
err_cnt += e;


np02_bmp_file_init_params bmp_init_params;
const double pixel_num = 8.0;
bmp_init_params.width_px = 320;
bmp_init_params.height_px = 320;
np02_bmp_file bmp_file(bmp_init_params);
const np02_xy xy_min(-20.0,-20.0);

bmp_debug_file_params bmp_file_p;
bmp_file_p.m_bmp_file = &bmp_file;
bmp_file_p.m_header = "shape_test_0.bmp";
bmp_file_p.m_footer.erase();
bmp_file_p.m_xy_min = xy_min; 
bmp_file_p.m_pixel_num = pixel_num;
bmp_file_p.m_basic_length = 40.0;
bmp_file_p.m_err_xy_found = false;
bmp_file_p.m_draw_overlap_circle = true;
bmp_file_p.m_draw_connecting_line = true;
bmp_file_p.m_shapes = &shape_vec;
bmp_file_p.m_color_vec = &color_vec;
bmp_file_p.m_err_xy_a = NULL;
bmp_file_p.m_err_xy_b = NULL;
bmp_file_p.m_err_shape_a = NULL;
bmp_file_p.m_err_shape_b = NULL;
bmp_file_p.m_debug_str_out = NULL;
draw_debug_shapes( bmp_file_p );

bmp_file.write_file( bmp_file_p.m_header.c_str() );

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
np02_shape_vec shape_vec;
std::vector<np02_bmp_color> color_vec;
np02_shp_alloc *shp_alloc = NULL;
advance_rand();
int ww = 3 + (get_rand_int() % 6);
CF01_HASH_CONSISTENCY_CHECK( cf01_obj_hash( 0, ww ) );
advance_rand();
int hh = 3 + (get_rand_int() % 6);
CF01_HASH_CONSISTENCY_CHECK( cf01_obj_hash( cf01_obj_hash( 0, ww ), hh ) );
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
    e = shp_alloc->verify_data(m_temp_err_msg,
        TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos);
    err_cnt += e;
    if( e > 0 ){
        print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 
        }
    }
make_rand_shapes( shp_alloc, ww, hh, basic_length, &shape_vec, &color_vec );
CF01_HASH_CONSISTENCY_CHECK( cf01_obj_hash( cf01_obj_hash( cf01_obj_hash( 
    0, ww ), hh ), np02_shape_vec_hash( shape_vec ) ) );
CF01_HASH_CONSISTENCY_CHECK( ( ( NULL != shp_alloc ) ? shp_alloc->hash() : 0 ) +
    cf01_obj_hash( cf01_obj_hash( cf01_obj_hash( 
    0, ww ), hh ), np02_shape_vec_hash( shape_vec ) ) );

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
const np02_shape *err_shape_a = NULL;
const np02_shape *err_shape_b = NULL;

for( i = 0; i < shape_vec.size(); ++i){
    assert( NULL != shape_vec.at(i) );
    e = shape_vec.at(i)->verify_data(m_temp_err_msg,
        TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos);
    err_cnt += e;
    if((e > 0) && !err_xy_found ){
        err_xy_found = true;
        should_write_file = true;
        shape_vec.at(i)->get_bb(&err_xy_a, &err_xy_b );
        err_shape_a = shape_vec.at(i);
        }
    if( e > 0 ){
        print_temp_buf_to_test_buf( __FILE__, __LINE__ );
        if(!err_xy_found ){
            err_xy_found = true;
            should_write_file = true;
            shape_vec.at(i)->get_bb(&err_xy_a, &err_xy_b );
            err_shape_a = shape_vec.at(i);
            }
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
            shape_i, shape_j, xy_a, xy_b, d, m_temp_err_msg,
            TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos);
        err_cnt += e;
        if( e > 0 ){
            print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 
            if( !err_xy_found ){
                err_xy_found = true;
                should_write_file = true;
                err_xy_a=xy_a;
                err_xy_b=xy_b;
                err_shape_a = shape_i;
                err_shape_b = shape_j;
                }
            }
        }
    }

if(should_write_file){
    np02_bmp_file_init_params bmp_init_params;
    static const double pixels_per_basic_length = 256.0;
    const double pixel_num = pixels_per_basic_length/basic_length;
    bmp_init_params.width_px =
        static_cast<int32_t>(floor(ww * pixels_per_basic_length));
    bmp_init_params.height_px =
        static_cast<int32_t>(floor(hh * pixels_per_basic_length));
    np02_bmp_file bmp_file(bmp_init_params);

    char bmp_file_name_buf[255];
    sprintf(&(bmp_file_name_buf[0]), "shape_test1_%i.bmp", m_iteration );
    std::string bmp_file_name( bmp_file_name_buf ); 
    const std::string prog_str = np02_test_main::get_prog_name() +
        std::string(" ") + np02_test_main::get_version_str();

    std::string debug_str_out;
    bmp_debug_file_params bmp_file_p;
    bmp_file_p.m_bmp_file = &bmp_file;
    bmp_file_p.m_header = bmp_file_name;
    bmp_file_p.m_footer = prog_str;
    bmp_file_p.m_xy_min = xy_min; 
    bmp_file_p.m_pixel_num = pixel_num;
    bmp_file_p.m_basic_length = basic_length;
    bmp_file_p.m_err_xy_found = err_xy_found;
    bmp_file_p.m_draw_overlap_circle = true;
    bmp_file_p.m_draw_connecting_line = true;
    bmp_file_p.m_shapes = &shape_vec;
    bmp_file_p.m_color_vec = &color_vec;
    bmp_file_p.m_err_xy_a = &err_xy_a;
    bmp_file_p.m_err_xy_b = &err_xy_b;
    bmp_file_p.m_err_shape_a = err_shape_a;
    bmp_file_p.m_err_shape_b = err_shape_b;
    bmp_file_p.m_debug_str_out = &debug_str_out;
    draw_debug_shapes( bmp_file_p );
    bmp_file.write_file( &(bmp_file_name[0]) );

    if(err_xy_found){
        np02_snprintf(m_temp_err_msg, TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos,
            "%s", debug_str_out.c_str()); 
        print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 

        char focus_bmp_file_name_buf[255];
        sprintf(&(focus_bmp_file_name_buf[0]), "shape_test1_%i_focus.bmp", m_iteration );
        std::string focus_bmp_file_name( focus_bmp_file_name_buf ); 

        np02_bmp_file_init_params focus_bmp_init_params;
        focus_bmp_init_params.width_px =
            static_cast<int32_t>(floor(ww * pixels_per_basic_length));
        focus_bmp_init_params.height_px =
            static_cast<int32_t>(floor(hh * pixels_per_basic_length));
        np02_bmp_file focus_bmp_file(focus_bmp_init_params);

        bmp_debug_file_params focus_bmp_file_p = bmp_file_p;
        focus_bmp_file_p.m_bmp_file = &focus_bmp_file;
        focus_bmp_file_p.m_header = focus_bmp_file_name;
        focus_bmp_file_p.m_shapes = NULL;
        focus_bmp_file_p.m_color_vec = NULL;
        draw_debug_shapes( focus_bmp_file_p );

        focus_bmp_file.write_file( &(focus_bmp_file_name[0]) );
        }
    }

CF01_HASH_CONSISTENCY_CHECK( ( ( NULL != shp_alloc ) ? shp_alloc->hash() : 0 ) +
    cf01_obj_hash( cf01_obj_hash( cf01_obj_hash( 
    0, ww ), hh ), np02_shape_vec_hash( shape_vec ) ) );

free_shapes( shp_alloc, &shape_vec );
if( NULL != shp_alloc ){ delete shp_alloc;  shp_alloc = NULL; }

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
const time_t start_time = time(NULL);
time_t prev_printf_time = 0;
for(; m_iteration < m_shape_test_iteration_count; ++m_iteration ){
    /* message */
    time_t now = time(NULL);
    const time_t printf_time_step = 120;
    if( now > ( prev_printf_time + printf_time_step ) ){
        const time_t elapsed_time = now - start_time;
        std::cout << "iteration:" << m_iteration
            << "   elapsed_time:" << elapsed_time << " s\n";
        prev_printf_time = now;
        }

    /* run test */
    switch(m_iteration % 2){
        default:
        case 0: err_cnt += execute_test_2_0(); break;
        case 1: err_cnt += execute_test_2_1(); break;
        }
    }
std::cout << "errors: " << err_cnt << "\n";
if( err_cnt > 0 ){ print_and_clear_test_buf(); }
if( NULL != AA_ERR_BUF_POS_PTR() && *(AA_ERR_BUF_POS_PTR()) > 0 ){ 
    std::cout << "\n\nAA_ERR_BUF:" << AA_ERR_BUF();
    }
AA_DECR_CALL_DEPTH();
return err_cnt;
}


int np02_shape_test::execute_test_2_0(){ 
AA_INCR_CALL_DEPTH();
int e = 0;
int err_cnt = 0;
int local_err_cnt = 0;
size_t i,j;
np02_shape_vec shape_vec;
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
    e = shp_alloc->verify_data(m_temp_err_msg,
        TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos);
    err_cnt += e;
    if( e > 0 ){
        print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 
        }
    }
make_rand_shapes( shp_alloc, ww, hh, basic_length, &shape_vec, &color_vec );
const np02_xy xy_min(0.0,0.0);

CF01_HASH_CONSISTENCY_CHECK( cf01_obj_hash( np02_shape_vec_hash( shape_vec ),
    ( NULL != shp_alloc ) ? shp_alloc->hash() : 0 ) );

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
    }
else if(50 > (get_rand_int() % 100)){
    loc_grid = new np02_loc_grid();
    }
else{
    loc_grid = &loc_grid_local;
    }

e = loc_grid->verify_data(m_temp_err_msg,
        TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos);
err_cnt += e;
if( e > 0 ){
    print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 
    }

loc_grid->set_extra_search_d(extra_search_d);
loc_grid->init_loc_grid(loc_grid_dim);

e = loc_grid->verify_data(m_temp_err_msg,
        TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos);
err_cnt += e;
if( e > 0 ){
    print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 
    }

CF01_HASH_CONSISTENCY_CHECK( cf01_obj_hash( loc_grid->hash(),
    cf01_obj_hash( np02_shape_vec_hash( shape_vec ),
    ( NULL != shp_alloc ) ? shp_alloc->hash() : 0 ) ) );

/* Put shapes in locator grid */
bool err_xy_found = false;
np02_xy err_xy_a(0.0, 0.0);
np02_xy err_xy_b(0.0, 0.0);
const np02_shape *err_shape_a = NULL;
const np02_shape *err_shape_b = NULL;
np02_shape *shape = NULL;
for( i = 0; i < shape_vec.size(); ++i){
    shape = shape_vec.at(i);
    e = shape->verify_data(m_temp_err_msg,
        TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos);
    err_cnt += e;
    if( e > 0 ){
        print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 
        if(!err_xy_found ){
            err_xy_found = true;
            should_write_file = true;
            shape->get_bb(&err_xy_a, &err_xy_b );
            err_shape_a = shape;
            err_shape_b = NULL;
            }
        }

    loc_grid->insert_shape_in_loc_grid(shape);
    
    e = shape->verify_data(m_temp_err_msg,
        TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos);
    err_cnt += e;
    if( e > 0 ){
        print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 
        if(!err_xy_found ){
            err_xy_found = true;
            should_write_file = true;
            shape->get_bb(&err_xy_a, &err_xy_b );
            err_shape_a = shape;
            }
        }

    CF01_HASH_CONSISTENCY_CHECK( cf01_obj_hash( loc_grid->hash(),
        cf01_obj_hash( np02_shape_vec_hash( shape_vec ),
        ( NULL != shp_alloc ) ? shp_alloc->hash() : 0 ) ) );
    }

e = loc_grid->verify_data(m_temp_err_msg,
        TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos);
AA_ALWAYS_ASSERT(0 == e);
err_cnt += e;
if( e > 0 ){
    print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 
    }

CF01_HASH_CONSISTENCY_CHECK( cf01_obj_hash( loc_grid->hash(),
    cf01_obj_hash( np02_shape_vec_hash( shape_vec ),
    ( NULL != shp_alloc ) ? shp_alloc->hash() : 0 ) ) );

/* 
From each shape,
  find nearest distance to every other shape
  find all local shapes with locator grid
  check that all shapes within locator grid extra-search-distance are found
    in the local locator grid search.
*/
np02_shape_vec local_shapes;
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
                np02_shape_vec_citr search_itr =
                    std::lower_bound( local_shapes.begin(), local_shapes.end(),
                        const_cast<np02_shape *>(shape_j));
                if( (search_itr == local_shapes.end()) || 
                    (shape_j != *search_itr)){
                    np02_snprintf(m_temp_err_msg, TEMP_ERR_MSG_CAP, 
                        &m_temp_err_msg_pos,
                        "shape_j not found in local_shapes"); 
                    print_temp_buf_to_test_buf( __FILE__, __LINE__ );
                    if(!err_xy_found ){
                        err_xy_found = true;
                        should_write_file = true;
                        shape_i->get_bb(&err_xy_a, &err_xy_b );
                        err_shape_a = shape_i;
                        err_shape_b = shape_j;
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
        static_cast<int32_t>(floor(ww * pixels_per_basic_length));
    bmp_init_params.height_px =
        static_cast<int32_t>(floor(hh * pixels_per_basic_length));
    np02_bmp_file bmp_file(bmp_init_params);

    loc_grid->write_bmp_file(xy_min, pixel_num, np02_bmp_color::gray(),
        &bmp_file);

    char bmp_file_name_buf[255];
    sprintf(&(bmp_file_name_buf[0]), "shape_test2_%i.bmp", m_iteration );
    std::string bmp_file_name( bmp_file_name_buf ); 
    const std::string prog_str = np02_test_main::get_prog_name() +
        std::string(" ") + np02_test_main::get_version_str();

    std::string debug_str_out;
    bmp_debug_file_params bmp_file_p;
    bmp_file_p.m_bmp_file = &bmp_file;
    bmp_file_p.m_header = bmp_file_name;
    bmp_file_p.m_footer = prog_str;
    bmp_file_p.m_xy_min = xy_min; 
    bmp_file_p.m_pixel_num = pixel_num;
    bmp_file_p.m_basic_length = basic_length;
    bmp_file_p.m_err_xy_found = err_xy_found;
    bmp_file_p.m_draw_overlap_circle = false;
    bmp_file_p.m_draw_connecting_line = false;
    bmp_file_p.m_shapes = &shape_vec;
    bmp_file_p.m_color_vec = &color_vec;
    bmp_file_p.m_err_xy_a = &err_xy_a;
    bmp_file_p.m_err_xy_b = &err_xy_b;
    bmp_file_p.m_err_shape_a = err_shape_a;
    bmp_file_p.m_err_shape_b = err_shape_b;
    bmp_file_p.m_debug_str_out = &debug_str_out;
    draw_debug_shapes( bmp_file_p );
    bmp_file.write_file( &(bmp_file_name[0]) );

    if(err_xy_found){
        np02_snprintf(m_temp_err_msg, TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos,
            "%s", debug_str_out.c_str()); 
        print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 

        char focus_bmp_file_name_buf[255];
        sprintf(&(focus_bmp_file_name_buf[0]), "shape_test2_%i_focus.bmp", m_iteration );
        std::string focus_bmp_file_name( focus_bmp_file_name_buf ); 

        np02_bmp_file_init_params focus_bmp_init_params;
        focus_bmp_init_params.width_px =
            static_cast<int32_t>(floor(ww * pixels_per_basic_length));
        focus_bmp_init_params.height_px =
            static_cast<int32_t>(floor(hh * pixels_per_basic_length));
        np02_bmp_file focus_bmp_file(focus_bmp_init_params);

        bmp_debug_file_params focus_bmp_file_p = bmp_file_p;
        focus_bmp_file_p.m_bmp_file = &focus_bmp_file;
        focus_bmp_file_p.m_header = focus_bmp_file_name;
        focus_bmp_file_p.m_shapes = NULL;
        focus_bmp_file_p.m_color_vec = NULL;
        draw_debug_shapes( focus_bmp_file_p );

        focus_bmp_file.write_file( &(focus_bmp_file_name[0]) );
        }
    }

err_xy_found = false;

/* 
Move some shapes
*/
for( i = 0; i < shape_vec.size(); ++i){
    np02_shape *shape_i = shape_vec.at(i);
    local_err_cnt = 0;
    e = shape_i->verify_data(m_temp_err_msg,
        TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos);
    AA_ALWAYS_ASSERT(0 == e);
    local_err_cnt += e;
    if( e > 0 ){
        print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 
        }

    if(shape_i->get_loc_grid() != loc_grid){
        ++local_err_cnt;
        np02_snprintf(m_temp_err_msg, TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos,
            "shape_i->get_loc_grid() != loc_grid" ); 
        print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 
        }

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

            e = shape_i->verify_data(m_temp_err_msg,
                TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos);
            local_err_cnt += e;
            if( e > 0 ){
                print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 
                }

            if(shape_i->get_loc_grid() != NULL){
                ++local_err_cnt;
                np02_snprintf(m_temp_err_msg, TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos,
                    "shape_i->get_loc_grid() != NULL" ); 
                print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 
                }
        

            shape_i->rotate_no_loc_grid(
                np02_xy(rot_ctr_x, rot_ctr_y), rot_deg);

            e = shape_i->verify_data(m_temp_err_msg,
                TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos);
            local_err_cnt += e;
            if( e > 0 ){
                print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 
                }

            loc_grid->insert_shape_in_loc_grid(shape_i);
            }
            break;
        case 9:
        case 10:
        case 11:
            {
            /* translate, explicitly update locator grid */
            loc_grid->remove_shape_from_loc_grid(shape_i);

            e = shape_i->verify_data(m_temp_err_msg,
                TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos);
            local_err_cnt += e;
            if( e > 0 ){
                print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 
                }

            if(shape_i->get_loc_grid() != NULL){
                ++local_err_cnt;
                np02_snprintf(m_temp_err_msg, TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos,
                    "shape_i->get_loc_grid() != NULL" ); 
                print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 
                }

            shape_i->translate_no_loc_grid(
                np02_xy(translation_dx, translation_dy));

            e = shape_i->verify_data(m_temp_err_msg,
                TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos);
            local_err_cnt += e;

            loc_grid->insert_shape_in_loc_grid(shape_i);
            }
            break;
        }

    e = shape_i->verify_data(m_temp_err_msg,
        TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos);
    local_err_cnt += e;
    if( e > 0 ){
        print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 
        }

    if(shape_i->get_loc_grid() != loc_grid){
        ++local_err_cnt;
        np02_snprintf(m_temp_err_msg, TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos,
            "shape_i->get_loc_grid() != loc_grid" ); 
        print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 
        }

    err_cnt += local_err_cnt;
    if((local_err_cnt > 0) && !err_xy_found ){
        err_xy_found = true;
        should_write_file = true;
        shape_i->get_bb(&err_xy_a, &err_xy_b );
        }

    CF01_HASH_CONSISTENCY_CHECK( cf01_obj_hash( loc_grid->hash(),
        cf01_obj_hash( np02_shape_vec_hash( shape_vec ),
        ( NULL != shp_alloc ) ? shp_alloc->hash() : 0 ) ) );
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
                np02_shape_vec_citr search_itr =
                    std::lower_bound( local_shapes.begin(), local_shapes.end(),
                        const_cast<np02_shape *>(shape_j));
                if( (search_itr == local_shapes.end()) || 
                    (shape_j != *search_itr)){
                    np02_snprintf(m_temp_err_msg, TEMP_ERR_MSG_CAP, 
                        &m_temp_err_msg_pos,
                        "shape_j not found in local_shapes"); 
                    print_temp_buf_to_test_buf( __FILE__, __LINE__ );
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
        static_cast<int32_t>(floor(ww * pixels_per_basic_length));
    bmp_init_params.height_px =
        static_cast<int32_t>(floor(hh * pixels_per_basic_length));
    np02_bmp_file bmp_file(bmp_init_params);

    loc_grid->write_bmp_file(xy_min, pixel_num, np02_bmp_color::gray(),
        &bmp_file);

    char bmp_file_name_buf[255];
    sprintf(&(bmp_file_name_buf[0]), "shape_test2a_%i.bmp", m_iteration );
    std::string bmp_file_name( bmp_file_name_buf ); 
    const std::string prog_str = np02_test_main::get_prog_name() +
        std::string(" ") + np02_test_main::get_version_str();

    std::string debug_str_out;
    bmp_debug_file_params bmp_file_p;
    bmp_file_p.m_bmp_file = &bmp_file;
    bmp_file_p.m_header = bmp_file_name;
    bmp_file_p.m_footer = prog_str;
    bmp_file_p.m_xy_min = xy_min; 
    bmp_file_p.m_pixel_num = pixel_num;
    bmp_file_p.m_basic_length = basic_length;
    bmp_file_p.m_err_xy_found = err_xy_found;
    bmp_file_p.m_draw_overlap_circle = false;
    bmp_file_p.m_draw_connecting_line = false;
    bmp_file_p.m_shapes = &shape_vec;
    bmp_file_p.m_color_vec = &color_vec;
    bmp_file_p.m_err_xy_a = &err_xy_a;
    bmp_file_p.m_err_xy_b = &err_xy_b;
    bmp_file_p.m_err_shape_a = err_shape_a;
    bmp_file_p.m_err_shape_b = err_shape_b;
    bmp_file_p.m_debug_str_out = &debug_str_out;
    draw_debug_shapes( bmp_file_p );
    bmp_file.write_file( &(bmp_file_name[0]) );

    if(err_xy_found){
        np02_snprintf(m_temp_err_msg, TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos,
            "%s", debug_str_out.c_str()); 
        print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 

        char focus_bmp_file_name_buf[255];
        sprintf(&(focus_bmp_file_name_buf[0]), "shape_test2a_%i_focus.bmp", m_iteration );
        std::string focus_bmp_file_name( focus_bmp_file_name_buf ); 

        np02_bmp_file_init_params focus_bmp_init_params;
        focus_bmp_init_params.width_px =
            static_cast<int32_t>(floor(ww * pixels_per_basic_length));
        focus_bmp_init_params.height_px =
            static_cast<int32_t>(floor(hh * pixels_per_basic_length));
        np02_bmp_file focus_bmp_file(focus_bmp_init_params);

        bmp_debug_file_params focus_bmp_file_p = bmp_file_p;
        focus_bmp_file_p.m_bmp_file = &focus_bmp_file;
        focus_bmp_file_p.m_header = focus_bmp_file_name;
        focus_bmp_file_p.m_shapes = NULL;
        focus_bmp_file_p.m_color_vec = NULL;
        draw_debug_shapes( focus_bmp_file_p );

        focus_bmp_file.write_file( &(focus_bmp_file_name[0]) );
        }
    }

e = loc_grid->verify_data(m_temp_err_msg,
    TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos);
    err_cnt += e;
if( e > 0 ){
    print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 
    }

if( NULL != shp_alloc ){
    e = shp_alloc->verify_data(m_temp_err_msg,
        TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos);
    err_cnt += e;
    if( e > 0 ){
        print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 
        }
    }
CF01_HASH_CONSISTENCY_CHECK( cf01_obj_hash( loc_grid->hash(),
    cf01_obj_hash( np02_shape_vec_hash( shape_vec ),
    ( NULL != shp_alloc ) ? shp_alloc->hash() : 0 ) ) );

free_shapes( shp_alloc, &shape_vec );

e = loc_grid->verify_data(m_temp_err_msg,
    TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos);
    err_cnt += e;
if( e > 0 ){
    print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 
    }

if( NULL != shp_alloc ){
    e = shp_alloc->verify_data(m_temp_err_msg,
        TEMP_ERR_MSG_CAP, &m_temp_err_msg_pos);
    err_cnt += e;
    if( e > 0 ){
        print_temp_buf_to_test_buf( __FILE__, __LINE__ ); 
        }
    }
CF01_HASH_CONSISTENCY_CHECK( cf01_obj_hash( loc_grid->hash(),
    cf01_obj_hash( np02_shape_vec_hash( shape_vec ),
    ( NULL != shp_alloc ) ? shp_alloc->hash() : 0 ) ) );

if( NULL != shp_alloc ){ delete shp_alloc; shp_alloc=NULL; }
loc_grid = NULL;

AA_DECR_CALL_DEPTH();

return err_cnt;
}

int np02_shape_test::execute_test_2_1(){ return execute_test_2_0(); }

int np02_shape_test::execute_test_3(){ return 0; }
int np02_shape_test::execute_test_4(){ return 0; }
int np02_shape_test::execute_test_5(){ return 0; }



} /* namespace np02 */
