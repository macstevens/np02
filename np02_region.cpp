/* np02_region.cpp  Newport 02 Library

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

np02_boundary_seg::np02_boundary_seg(): np02_shape_owner(), m_shp_alloc(NULL),
    m_alloc_idx(0), m_owner(NULL), m_shape(NULL), m_prev(NULL),
    m_next(NULL), m_orientation(NP02_BOUNDARY_SEG_ORIENTATION_COUNT),
    m_invert_status(NP02_BOUNDARY_SEG_INVERT_STATUS_COUNT){
}

np02_boundary_seg::~np02_boundary_seg(){
destruct();
}

void np02_boundary_seg::destruct(){
clear_shape();
}

void np02_boundary_seg::set_shape(np02_shape *shape){
if( m_shape != shape ){
    clear_shape();
    }
m_shape = shape;
}

np02_shape *np02_boundary_seg::get_shape() const{
return m_shape;
}

uint64_t np02_boundary_seg::hash( const uint64_t& h_in ) const{
/*TODO: implement */
return 0;
}

int np02_boundary_seg::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const{
int err_cnt = 0;
err_cnt += np02_shape_owner::verify_data(err_msg,err_msg_capacity,err_msg_pos);
/* TODO: check m_shp_alloc contains this */
/* TODO: check m_shp_alloc matches m_shape */
/* TODO: check m_shape points to this */

if( NULL != m_prev ){
    if( this != m_prev->m_next ){
        ++err_cnt;
        /* TODO: message*/
        }
    if( m_owner != m_prev->m_owner ){
        ++err_cnt;
        /* TODO: message*/
        }
    }

if( NULL != m_next ){
    if( this != m_next->m_prev ){
        ++err_cnt;
        /* TODO: message*/
        }
    if( m_owner != m_next->m_owner ){
        ++err_cnt;
        /* TODO: message*/
        }
    }

/* TODO: check m_orientation matches shape type */
/* TODO: check m_invert_status matches shape type */
return err_cnt;
}

void np02_boundary_seg::clear_shape(){
if( NULL != m_shape ){
    AUTO_ASSERT( m_shape->get_shp_alloc() == m_shp_alloc );
    if(NULL == m_shp_alloc){
        delete m_shape;
        }
    else{
        m_shp_alloc->free_shape(m_shape);
        }
    m_shape = NULL;
    }
}


np02_boundary::np02_boundary(): m_shp_alloc(NULL), m_alloc_idx(0),
    m_owner(NULL), m_segs_head(NULL), m_segs_tail(NULL){
}

np02_boundary::~np02_boundary(){
destruct();
}

void np02_boundary::destruct(){
clear_boundary_segs();
}

void np02_boundary::clear_boundary_segs(){
/* TODO: implement */
}

void np02_boundary::add_boundary_seg(np02_boundary_seg *s){
/* TODO: implement */
}

uint64_t np02_boundary::hash( const uint64_t& h_in ) const{
/*TODO: implement */
return 0;
}

int np02_boundary::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const{
int err_cnt = 0;
return err_cnt;
}


np02_region::np02_region(): m_shp_alloc(NULL), m_alloc_idx(0),
    m_boundaries(){
}

np02_region::~np02_region(){
destruct();
}

void np02_region::destruct(){
clear_boundaries();
}

void np02_region::clear_boundaries(){
const size_t boundary_count = m_boundaries.size();
for( size_t i = 0; i < boundary_count; ++i ){
    np02_boundary *b = m_boundaries.at(i);
    AA_ALWAYS_ASSERT(NULL != b)
    AUTO_ASSERT( b->get_shp_alloc() == m_shp_alloc );
    AUTO_ASSERT( b->get_owner() == this );
    if( NULL == m_shp_alloc ){
        delete b;
        }
    else{
        /* TODO: m_shp_alloc->free_boundary(b); */
        }
    }
}

void np02_region::add_boundary(np02_boundary *b){
AA_ALWAYS_ASSERT(NULL != b)
AUTO_ASSERT( b->get_shp_alloc() == m_shp_alloc );
m_boundaries.push_back(b);
b->set_owner(this);
}

uint64_t np02_region::hash( const uint64_t& h_in ) const{
/*TODO: implement */
return 0;
}

int np02_region::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const{
int err_cnt = 0;
/* TODO: boundaries non-NULL */
/* TODO: boundaries unique */
/* TODO: this == boundary->owner */
/* TODO: boundary->verify_data() */


return err_cnt;
}




} /* namespace np02 */
