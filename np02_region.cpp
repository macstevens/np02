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
    #define AA_CALL_DEPTH_BLOCK()       CF01_AA_CALL_DEPTH_BLOCK()
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
    #define AA_CALL_DEPTH_BLOCK()
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
    m_next(NULL){
}

np02_boundary_seg::~np02_boundary_seg(){
destruct();
}

void np02_boundary_seg::destruct(){
clear_shape();
}

np02_shp_alloc *np02_boundary_seg::get_shp_alloc() const{
return m_shp_alloc;
}

np02_region *np02_boundary_seg::get_region() const{
np02_region *region = (NULL == m_owner) ? NULL : m_owner->get_owner();
return region;
}

void np02_boundary_seg::set_shape(np02_shape *shape){
if( m_shape != shape ){
    clear_shape();
    m_shape = shape;
    if( NULL != m_shape){
        m_shape->set_shape_owner(this);
        }
    }
AUTO_ASSERT((NULL == m_shape) || (this == m_shape->get_boundary_seg()));
}

np02_shape *np02_boundary_seg::get_shape() const{
AUTO_ASSERT((NULL == m_shape) || (this == m_shape->get_boundary_seg()));
return m_shape;
}

np02_loc_grid *np02_boundary_seg::get_loc_grid() const{
np02_loc_grid *loc_grid = 
    (NULL == m_shape) ? NULL : m_shape->get_loc_grid();
return loc_grid;
}

void np02_boundary_seg::insert_boundary_seg_in_loc_grid(
    np02_loc_grid *loc_grid){
remove_boundary_seg_from_loc_grid();
if( ( NULL != loc_grid ) && ( NULL != m_shape ) ){
    loc_grid->insert_shape_in_loc_grid(m_shape);
    }
AUTO_ASSERT( ( NULL == m_shape ) || ( get_loc_grid() == loc_grid ) );
}

void np02_boundary_seg::remove_boundary_seg_from_loc_grid(){
np02_loc_grid *loc_grid = get_loc_grid();
if( ( NULL != loc_grid ) && ( NULL != m_shape ) ){
    loc_grid->remove_shape_from_loc_grid(m_shape);
    }
AUTO_ASSERT( NULL == get_loc_grid() );
}

np02_xy np02_boundary_seg::get_p0() const{
np02_xy p0(0.0,0.0);
if(NULL != m_shape){
    const np02_shape_orientation orientation=m_shape->get_shape_orientation();
    const np02_circle *c = dynamic_cast<const np02_circle *>(m_shape);
    const np02_arc *a = dynamic_cast<const np02_arc *>(m_shape);
    const np02_line_seg *n = dynamic_cast<const np02_line_seg *>(m_shape);
    const np02_rect *r = dynamic_cast<const np02_rect *>(m_shape);
    const np02_polygon *p = dynamic_cast<const np02_polygon *>(m_shape);
    const np02_spline *i = dynamic_cast<const np02_spline *>(m_shape);
    if(NULL != c){}
    else if(NULL!=a){
        p0 = (NP02_SHAPE_ORIENTATION_FWD == orientation )?
            a->get_p_0() : a->get_p_1();
        }
    else if(NULL!=n){
        p0 = (NP02_SHAPE_ORIENTATION_FWD == orientation )?
            n->get_p_0() : n->get_p_1();
        }
    else if(NULL!=r){}
    else if(NULL!=p){ }
    else if(NULL!=i){ AA_ALWAYS_ASSERT(false); /* spline query not supported */ }
    else{ AA_ALWAYS_ASSERT(false); }
    }
return p0;
}

np02_xy np02_boundary_seg::get_p1() const{
np02_xy p1(0.0,0.0);
if(NULL != m_shape){
    const np02_shape_orientation orientation=m_shape->get_shape_orientation();
    const np02_circle *c = dynamic_cast<const np02_circle *>(m_shape);
    const np02_arc *a = dynamic_cast<const np02_arc *>(m_shape);
    const np02_line_seg *n = dynamic_cast<const np02_line_seg *>(m_shape);
    const np02_rect *r = dynamic_cast<const np02_rect *>(m_shape);
    const np02_polygon *p = dynamic_cast<const np02_polygon *>(m_shape);
    const np02_spline *i = dynamic_cast<const np02_spline *>(m_shape);
    if(NULL != c){}
    else if(NULL!=a){
        p1 = (NP02_SHAPE_ORIENTATION_FWD == orientation )?
            a->get_p_1() : a->get_p_0();
        }
    else if(NULL!=n){
        p1 = (NP02_SHAPE_ORIENTATION_FWD == orientation )?
            n->get_p_1() : n->get_p_0();
        }
    else if(NULL!=r){}
    else if(NULL!=p){ }
    else if(NULL!=i){ AA_ALWAYS_ASSERT(false); /* spline query not supported */ }
    else{ AA_ALWAYS_ASSERT(false); }
    }
return p1;
}

double np02_boundary_seg::get_boundary_gap_to_next() const{
double gap = 0.0;
if( NULL != m_next ){
    const np02_xy this_p1 = get_p1();
    const np02_xy next_p0 = m_next->get_p0();
    gap = this_p1.get_distance_to(next_p0);    
    }
return gap;
}

void np02_boundary_seg::translate(const np02_xy& dxy){
if( NULL != m_shape ){
    m_shape->translate(dxy);
    }
}

void np02_boundary_seg::translate_no_loc_grid(const np02_xy& dxy){
if( NULL != m_shape ){
    m_shape->translate_no_loc_grid(dxy);
    }
}

void np02_boundary_seg::rotate(const np02_xy& rot_ctr, const double& rot_deg){
if( NULL != m_shape ){
    m_shape->rotate(rot_ctr, rot_deg);
    }
}

void np02_boundary_seg::rotate_no_loc_grid(const np02_xy& rot_ctr,
    const double& rot_deg){
if( NULL != m_shape ){
    m_shape->rotate_no_loc_grid(rot_ctr, rot_deg);
    }
}

void np02_boundary_seg::invert(){
AA_CALL_DEPTH_BLOCK();
if( NULL != m_shape ){
    /* TODO: make np02_shape::invert() function  
     (maybe eliminate np02_shape_orientation flag because it
     is equivalent to np02_shape_invert_status)*/
    np02_circle *c = dynamic_cast<np02_circle *>(m_shape);
    np02_arc *a = dynamic_cast<np02_arc *>(m_shape);
    np02_line_seg *n = dynamic_cast<np02_line_seg *>(m_shape);
    np02_rect *r = dynamic_cast<np02_rect *>(m_shape);
    np02_polygon *p = dynamic_cast<np02_polygon *>(m_shape);
    np02_spline *i = dynamic_cast<np02_spline *>(m_shape);
    if(NULL != c){ c->invert(); }
    else if(NULL!=a){ a->reverse_orientation(); }
    else if(NULL!=n){ n->reverse_orientation(); }
    else if(NULL!=r){ r->invert(); }
    else if(NULL!=p){ p->invert(); }
    else if(NULL!=i){ AA_ALWAYS_ASSERT(false); /* spline query not supported */ }
    else{ AA_ALWAYS_ASSERT(false); }
    }
}

np02_boundary_seg *np02_boundary_seg::copy_boundary_seg(
    np02_shp_alloc *shp_alloc) const{
np02_boundary_seg *s = (NULL == shp_alloc) ?
    new np02_boundary_seg() : shp_alloc->alloc_boundary_seg();
s->init_from_boundary_seg(this);
return s;
}

void np02_boundary_seg::init_from_boundary_seg(
    const np02_boundary_seg *master){
AA_CALL_DEPTH_BLOCK();
if( (NULL != master) && (this != master) ){
    destruct();
    AUTO_ASSERT(NULL == m_owner);
    AUTO_ASSERT(NULL == m_shape);
    AUTO_ASSERT(NULL == m_prev);
    AUTO_ASSERT(NULL == m_next);
    const np02_shape *master_shape = master->m_shape;
    if( NULL != master_shape ){
        m_shape = master_shape->copy_shape(m_shp_alloc);
        m_shape->set_shape_owner(this);
        }
    AUTO_ASSERT( 0 == compare_boundary_seg(master) );
    }
}

int np02_boundary_seg::compare_boundary_seg(
    const np02_boundary_seg *other) const{
int result = 0;
if( NULL == other ){
    result = 1;
    }
else{
    const np02_shape *other_shape = other->m_shape;
    if( NULL == m_shape ){
        if( NULL == other_shape ){
            AUTO_ASSERT(0 == result);
            }
        else{
            result = -1;
            }
        }
    else if( NULL == other_shape ){
        result = 1;
        }
    else{
        result = m_shape->compare(other_shape);
        }
    }
return result;
}

uint64_t np02_boundary_seg::hash( const uint64_t& h_in ) const{
uint64_t h = h_in;
#if defined( CF01_SUPPORT )
h = cf01_obj_hash( h, m_alloc_idx );
#else
h += reinterpret_cast<cf01_uint64>(m_alloc_idx);
h ^= ((h << 27) | (h >> 37));
#endif
if( NULL != m_shape ){
    h = m_shape->hash(h);
    }
return h;
}

double np02_boundary_seg::get_distance_from_xy(const np02_xy& xy,
    np02_xy *near_xy, double *d_from2, bool *valid )const{
double d_from = 0.0;
bool vld = false;
if( NULL != m_shape ){
    d_from = m_shape->get_bseg_dist_from_xy( xy, near_xy, d_from2 );
    vld = true;
    }
if( NULL != valid ){
    *valid = vld;
    }
return d_from;
}

int np02_boundary_seg::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const{
int err_cnt = 0;
err_cnt += np02_shape_owner::verify_data(err_msg,err_msg_capacity,err_msg_pos);
if( NULL != m_shp_alloc ){
    const np02_boundary_seg *b_seg =
        m_shp_alloc->alloc_get_boundary_seg_by_idx(m_alloc_idx);
    if( this != b_seg ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary_seg: this=%llx  != "
            "(shp_alloc=%llx) boundary_seg_by_idx(idx=%i)=%llx\n",
            this, m_shp_alloc, m_alloc_idx, b_seg ); 
        }
    }

if( NULL != m_shape ){
    err_cnt += m_shape->verify_data(err_msg,err_msg_capacity,err_msg_pos);
    const np02_shape_owner *shp_ownr = m_shape->get_shape_owner();
    if( this != shp_ownr ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary_seg: this=%llx  != (shape=%llx) owner=%llx\n",
            this, m_shape, shp_ownr ); 
        }
    const np02_boundary_seg* shp_bdry_seg = m_shape->get_boundary_seg();
    if (this != shp_bdry_seg) {
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary_seg: this=%llx  != (shape=%llx) boundary_seg=%llx\n",
            this, m_shape, shp_bdry_seg);
        }
    const np02_shp_alloc *sa = m_shape->get_shp_alloc();
    if( m_shp_alloc != sa ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary_seg: (this=%llx) shp_alloc=%llx != "
            "(shape=%llx) shp_alloc=%llx\n",
            this, m_shp_alloc, m_shape, sa );
        }
    }

if( NULL != m_prev ){
    if( this != m_prev->m_next ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary_seg: this=%llx != (prev=%llx)->(next=%llx)\n",
            this, m_prev, m_prev->m_next ); 
        }
    if( m_owner != m_prev->m_owner ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary_seg: (this=%llx) (owner=%llx) != (prev=%llx)->(owner=%llx)\n",
            this, m_owner, m_prev, m_prev->m_owner ); 
        }
    }

if( NULL != m_next ){
    if( this != m_next->m_prev ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary_seg: this=%llx != (next=%llx)->(prev=%llx)\n",
            this, m_next, m_next->m_prev ); 
        }
    if( m_owner != m_next->m_owner ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary_seg: (this=%llx) (owner=%llx) != (next=%llx)->(owner=%llx)\n",
            this, m_owner, m_next, m_next->m_owner ); 
        }
    }
return err_cnt;
}

int np02_boundary_seg::verify_loc_grid_properly_contains_boundary_seg(
    const np02_loc_grid* loc_grid, const np02_boundary_seg* boundary_seg,
    char* err_msg, const size_t err_msg_capacity,
    size_t* err_msg_pos) {
int err_cnt = 0;
const np02_shape *shape = (NULL == boundary_seg) ? NULL :
    boundary_seg->get_shape();
if( NULL == boundary_seg ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "boundary_seg NULL\n" );
    }
else if (NULL == shape) {
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "boundary_seg=%llx : shape NULL\n", boundary_seg);
    }
else{
    err_cnt += np02_shape::verify_loc_grid_properly_contains_shape( loc_grid, 
        shape, err_msg, err_msg_capacity, err_msg_pos );
    }
return err_cnt;
}

std::ostream& np02_boundary_seg::ostream_output(std::ostream& os) const{
/* TODO: implement */
return os;
}

void np02_boundary_seg::write_bmp_file(const np02_xy& xy_min,
    const double& pixel_num, const np02_bmp_color& color,
    np02_bmp_file *bmp_file) const{
if( NULL != m_shape ){
    m_shape->write_bmp_file(xy_min, pixel_num, color, bmp_file);
    }
}

void np02_boundary_seg::write_bmp_file_edge(const np02_xy& xy_min,
    const double& pixel_num, const np02_bmp_color& color,
    np02_bmp_file *bmp_file) const{
if( NULL != m_shape ){
    m_shape->write_bmp_file_edge(xy_min, pixel_num, color, bmp_file);
    }
}

void np02_boundary_seg::write_dxf_file(const std::string& layer,
    const uint8_t& color, np02_dxf_file *dxf_file) const{
if( NULL != m_shape ){
    m_shape->write_dxf_file(layer, color, dxf_file);
    }
}

void np02_boundary_seg::clear_shape(){
if( NULL != m_shape ){
    AUTO_ASSERT( m_shape->get_shape_owner() == this );
    AUTO_ASSERT( m_shape->get_boundary_seg() == this );
    AUTO_ASSERT( m_shape->get_shp_alloc() == m_shp_alloc );
    m_shape->set_shape_owner(NULL);
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
    m_owner(NULL), m_segs_head(NULL), m_segs_tail(NULL), 
    m_rgn_bdry_prev(NULL), m_rgn_bdry_next(NULL){
}

np02_boundary::~np02_boundary(){
destruct();
}

void np02_boundary::destruct(){
if( NULL != m_owner ){
    m_owner->remove_boundary(this);
    }
clear_boundary_segs();
}

void np02_boundary::clear_boundary_segs(){
AUTO_ASSERT( 0 == verify_data_segs_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ))
while( NULL != m_segs_tail ){
    np02_boundary_seg *s = m_segs_tail;
    remove_boundary_seg(s);
    if( NULL == m_shp_alloc ){
        delete s;
        }
    else{
        m_shp_alloc->free_boundary_seg(s);
        }
    }
AUTO_ASSERT( 0 == verify_data_segs_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ))
}

/* add boundary_seg, make loop sequence */
void np02_boundary::add_boundary_seg_loop(np02_boundary_seg *s){
AA_CALL_DEPTH_BLOCK();
AUTO_ASSERT( 0 == verify_data_segs_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ))
np02_boundary *ownr = s->get_owner();
if( NULL != ownr ){
    ownr->remove_boundary_seg(s);
    }
AUTO_ASSERT( NULL == s->get_owner());
AUTO_ASSERT( NULL == s->get_prev());
AUTO_ASSERT( NULL == s->get_next());
s->set_owner(this);
if( NULL == m_segs_head){
    AUTO_ASSERT(NULL == m_segs_tail);
    m_segs_head = s;
    m_segs_tail = s;
    s->set_prev(s);
    s->set_next(s);
    }
else{
    AUTO_ASSERT(NULL != m_segs_tail);

    /* create loop, if chain is currently open */
    if(NULL == m_segs_tail->get_next()){
        AUTO_ASSERT(NULL == m_segs_head->get_prev());
        m_segs_head->set_prev(m_segs_tail);
        m_segs_tail->set_next(m_segs_head);
        m_segs_tail = m_segs_head;
        }
    AUTO_ASSERT(m_segs_head == m_segs_tail);
    AUTO_ASSERT(NULL != m_segs_head->get_prev());
    AUTO_ASSERT(NULL != m_segs_tail->get_next());

    /* insert before tail(head) */
    np02_boundary_seg *p = m_segs_tail->get_prev();
    p->set_next(s);
    s->set_prev(p);
    m_segs_tail->set_prev(s);
    s->set_next(m_segs_tail);
    }
AUTO_ASSERT(NULL != m_segs_head->get_prev());
AUTO_ASSERT(NULL != m_segs_tail->get_next());
AUTO_ASSERT(m_segs_head == m_segs_tail);
AUTO_ASSERT( 0 == verify_data_segs_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ));
AUTO_ASSERT( 0 == verify_data_segs_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ));
}

/* add boundary_seg, make open (non-loop) sequence */
void np02_boundary::add_boundary_seg_open(np02_boundary_seg *s){
AA_CALL_DEPTH_BLOCK();
AUTO_ASSERT( 0 == verify_data_segs_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ))
np02_boundary *ownr = s->get_owner();
if( NULL != ownr ){
    ownr->remove_boundary_seg(s);
    }
AUTO_ASSERT( NULL == s->get_owner());
AUTO_ASSERT( NULL == s->get_prev());
AUTO_ASSERT( NULL == s->get_next());
s->set_owner(this);
if( NULL == m_segs_head){
    AUTO_ASSERT(NULL == m_segs_tail);
    m_segs_head = s;
    m_segs_tail = s;
    }
else{
    AUTO_ASSERT(NULL != m_segs_tail);

    /* break loop, if chain is currently a loop */
    if(NULL != m_segs_tail->get_next()){
        AUTO_ASSERT(m_segs_head == m_segs_tail);
        m_segs_tail = m_segs_head->get_prev();
        m_segs_head->set_prev(NULL);
        m_segs_tail->set_next(NULL);
        }
    AUTO_ASSERT(NULL == m_segs_head->get_prev());
    AUTO_ASSERT(NULL == m_segs_tail->get_next());

    /* add to end */
    m_segs_tail->set_next(s);
    s->set_prev(m_segs_tail);
    m_segs_tail = s;
    }
AUTO_ASSERT(NULL == m_segs_head->get_prev());
AUTO_ASSERT(NULL == m_segs_tail->get_next());
AUTO_ASSERT(s == m_segs_tail);
AUTO_ASSERT(0 == verify_data_segs_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR()));
}

void np02_boundary::remove_boundary_seg(np02_boundary_seg *s){
AA_CALL_DEPTH_BLOCK();
AUTO_ASSERT( 0 == verify_data_segs_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ));
AUTO_ASSERT( (NULL != s ) && ( this == s->get_owner() ) );
if( (NULL != s ) && ( this == s->get_owner() ) ){
    AA_ALWAYS_ASSERT(NULL != m_segs_head);
    AA_ALWAYS_ASSERT(NULL != m_segs_tail);
    np02_boundary_seg *p = s->get_prev();
    np02_boundary_seg *n = s->get_next();
    if( p == s ){
        /* single element loop */
        AUTO_ASSERT(s == n);
        AUTO_ASSERT(s == m_segs_head);
        AUTO_ASSERT(s == m_segs_tail);
        m_segs_head = NULL;
        m_segs_tail = NULL;
        }
    else if( (m_segs_head == m_segs_tail) && (NULL == p) ){
        /* single element open list */
        AUTO_ASSERT(NULL == n);
        AUTO_ASSERT(s == m_segs_head);
        AUTO_ASSERT(s == m_segs_tail);
        m_segs_head = NULL;
        m_segs_tail = NULL;
        }
    else{
        if( m_segs_head == s ){
            AUTO_ASSERT(NULL != n);
            m_segs_head = n;
            if( m_segs_tail == s ){ 
                /* loop */
                m_segs_tail = n;
                }
            }
        else if( m_segs_tail == s ){ 
            /* open */
            AUTO_ASSERT(m_segs_head != m_segs_tail);
            AUTO_ASSERT(NULL != p);
            AUTO_ASSERT(NULL == m_segs_head->get_prev());
            AUTO_ASSERT(NULL != m_segs_head->get_next());
            AUTO_ASSERT(NULL != m_segs_tail->get_prev());
            AUTO_ASSERT(NULL == m_segs_tail->get_next());
            m_segs_tail = p;
            }
        if( NULL != p ){ p->set_next(n); }
        if( NULL != n ){ n->set_prev(p); }
        }
    s->set_owner(NULL);
    s->set_prev(NULL);
    s->set_next(NULL);
    }

AUTO_ASSERT( 0 == verify_data_segs_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ));
}

size_t np02_boundary::get_boundary_seg_count() const{
size_t c = 0;
const np02_boundary_seg *s = m_segs_head;
while( NULL != s ){
    ++c;
    s = advance_boundary_seg_citr(s);
    }
return c;
}

size_t np02_boundary::get_shape_count() const{
size_t c = 0;
const np02_boundary_seg *s = m_segs_head;
while( NULL != s ){
    if(NULL != s->get_shape()){
        ++c;
        }
    s = advance_boundary_seg_citr(s);
    }
return c;
}

/* stop counting when shape count reaches threshold */
void np02_boundary::get_shape_count_lower_bound(
    const size_t& max_shape_count_l_b, size_t *count ) const{
if( NULL != count ){
    const np02_boundary_seg *s = m_segs_head;
    while( ( NULL != s ) && ( *count < max_shape_count_l_b ) ) {
        if(NULL != s->get_shape()){
            ++*count;
            }
        s = advance_boundary_seg_citr(s);
        }
    }
}

np02_shape *np02_boundary::get_first_shape() const{
np02_shape *shape = NULL;
const np02_boundary_seg *s = m_segs_head;
while( (NULL != s) && (NULL == shape) ){
    shape = s->get_shape();
    s = advance_boundary_seg_citr(s);
    }
return shape;
}

np02_loc_grid *np02_boundary::get_loc_grid() const{
np02_shape *shape = get_first_shape();
np02_loc_grid *loc_grid = ( NULL == shape ) ? NULL : shape->get_loc_grid();
return loc_grid;
}

void np02_boundary::insert_boundary_in_loc_grid(np02_loc_grid *loc_grid){
np02_boundary_seg *s = m_segs_head;
while( NULL != s ){
    s->insert_boundary_seg_in_loc_grid(loc_grid);
    s = advance_boundary_seg_itr(s);
    }
AUTO_ASSERT( ( NULL == get_first_shape() ) || ( get_loc_grid() == loc_grid ) );
}

void np02_boundary::remove_boundary_from_loc_grid(){
np02_boundary_seg *s = m_segs_head;
while( NULL != s ){
    s->remove_boundary_seg_from_loc_grid();
    s = advance_boundary_seg_itr(s);
    }
AUTO_ASSERT( NULL == get_loc_grid() );
}

/* advance boundary seg pointer as an iterator.  returns NULL when iterator
   reaches end of list */
np02_boundary_seg *np02_boundary::advance_boundary_seg_itr(
    np02_boundary_seg *bndry_seg_itr) const{
AUTO_ASSERT( 0 == verify_data_segs_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ))
AUTO_ASSERT( (NULL == bndry_seg_itr) || (this == bndry_seg_itr->get_owner()) );
np02_boundary_seg *s = bndry_seg_itr;
if( NULL != s ){
    s = s->get_next();
    if( ( NULL != s ) && ( ( s == m_segs_head ) &&
        ( m_segs_head == m_segs_tail ) ) ){
        s = NULL; 
        }
    }
return s;
}


void np02_boundary::get_bb(np02_xy *xy_min, np02_xy *xy_max, bool *valid)const{
np02_xy xy_low(0.0,0.0);
np02_xy xy_high(0.0,0.0);
np02_xy s_xy_low(0.0,0.0);
np02_xy s_xy_high(0.0,0.0);
bool ok = false;
const np02_boundary_seg *s = m_segs_head;
while( NULL != s ){
    const np02_shape *shape = s->get_shape();
    if( NULL != shape ){
        if(ok){
            shape->get_bb(&s_xy_low,&s_xy_high);
            np02_shape::bb_merge(xy_low,xy_high,s_xy_low,s_xy_high,
                &xy_low,&xy_high );
            }
        else{
            shape->get_bb(&xy_low,&xy_high);
            ok = true;
            }
        }
    s = advance_boundary_seg_citr(s);
    }
if( NULL != xy_min ){
    *xy_min = xy_low;
    }
if( NULL != xy_max ){
    *xy_max = xy_high;
    }
if( NULL != valid ){
    *valid = ok;
    }
}

double np02_boundary::get_boundary_gap_sum() const{
double gap_sum = 0.0;
const np02_boundary_seg *s = m_segs_head;
if( ( NULL != m_segs_head ) && ( m_segs_head == m_segs_tail ) &&
    ( NULL != m_segs_head->get_prev() ) ){
    /* loop*/ 
    while( NULL != s ){
        gap_sum += s->get_boundary_gap_to_next();
        s = advance_boundary_seg_citr(s);
        }
    }
return gap_sum;
}

void np02_boundary::translate(const np02_xy& dxy){
np02_boundary_seg *s = m_segs_head;
while( NULL != s ){
    s->translate(dxy);
    s = advance_boundary_seg_itr(s);
    }
}

void np02_boundary::translate_no_loc_grid(const np02_xy& dxy){
np02_boundary_seg *s = m_segs_head;
while( NULL != s ){
    s->translate_no_loc_grid(dxy);
    s = advance_boundary_seg_itr(s);
    }
}

void np02_boundary::rotate(const np02_xy& rot_ctr, const double& rot_deg){
np02_boundary_seg *s = m_segs_head;
while( NULL != s ){
    s->rotate(rot_ctr, rot_deg);
    s = advance_boundary_seg_itr(s);
    }
}

void np02_boundary::rotate_no_loc_grid(const np02_xy& rot_ctr,
    const double& rot_deg){
np02_boundary_seg *s = m_segs_head;
while( NULL != s ){
    s->rotate_no_loc_grid(rot_ctr, rot_deg);
    s = advance_boundary_seg_itr(s);
    }
}

void np02_boundary::invert(){
if( NULL != m_segs_head ){
    AUTO_ASSERT( 0 == verify_data(AA_ERR_BUF(),
        AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ));

    /* invert each boundary segment */
    np02_boundary_seg *s = m_segs_head;
    while( NULL != s ){
        s->invert();
        s = advance_boundary_seg_itr(s);
        }
    
    /* reverse sequence */
    s = m_segs_head;
    while( NULL != s ){
        np02_boundary_seg *p = s->get_prev();
        np02_boundary_seg *n = s->get_next();
        s->set_prev(n);
        s->set_next(p);
        s = n;
        if( ( NULL != s ) && ( ( s == m_segs_head ) &&
            ( m_segs_head == m_segs_tail ) ) ){
            s = NULL; 
            }        
        }

    /* swap head & tail */
    if( m_segs_head != m_segs_tail){ 
        s = m_segs_head;
        m_segs_head = m_segs_tail;
        m_segs_tail = s;
        }

    AUTO_ASSERT( 0 == verify_data(AA_ERR_BUF(),
        AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ));
    }
}

np02_boundary *np02_boundary::copy_boundary(np02_shp_alloc *shp_alloc) const{
np02_boundary *b = (NULL == shp_alloc) ?
    new np02_boundary() : shp_alloc->alloc_boundary();
b->init_from_boundary(this);
return b;
}

void np02_boundary::init_from_boundary( const np02_boundary *master ){
AA_CALL_DEPTH_BLOCK();
if( (NULL != master) && (this != master) ){
    AUTO_ASSERT(0 == verify_data(AA_ERR_BUF(), AA_ERR_BUF_CAPACITY(),
        AA_ERR_BUF_POS_PTR()));
    AA_ALWAYS_ASSERT(0 == master->verify_data(AA_ERR_BUF(), AA_ERR_BUF_CAPACITY(),
        AA_ERR_BUF_POS_PTR()));
    AUTO_ASSERT(0 == master->verify_data(AA_ERR_BUF(), AA_ERR_BUF_CAPACITY(),
        AA_ERR_BUF_POS_PTR()));
    clear_boundary_segs();
    const np02_boundary_seg *s = master->m_segs_head;
    const bool is_loop = (( NULL != s ) && ( NULL != s->get_prev() )) ?
        true : false;
    AUTO_ASSERT( (!is_loop) || ( master->m_segs_tail == s) );
    while( NULL != s ){
        np02_boundary_seg *s_copy = s->copy_boundary_seg(m_shp_alloc);
        if( is_loop ){ add_boundary_seg_loop(s_copy); }
        else{ add_boundary_seg_open(s_copy); }
        s = master->advance_boundary_seg_citr(s);
        }
    AUTO_ASSERT(0 == verify_data(AA_ERR_BUF(), AA_ERR_BUF_CAPACITY(),
        AA_ERR_BUF_POS_PTR()));
    AUTO_ASSERT( 0 == compare_boundary(master) );
    }
}

int np02_boundary::compare_boundary( const np02_boundary *other ) const{
int result = 0;
if( NULL == other ){
    result = 1;
    }
else{
    const np02_boundary_seg *s = m_segs_head;
    const np02_boundary_seg *other_s = other->m_segs_head;
    while( ( NULL != s ) && ( NULL != other_s ) && ( 0 == result ) ){
        result = s->compare_boundary_seg(other_s);
        s = advance_boundary_seg_citr(s);
        other_s = other->advance_boundary_seg_citr(other_s);
        }
    if( ( 0 == result ) && ( ( ( NULL != s ) || ( NULL != other_s ) ) ) ){
        if( NULL == s ){
            AUTO_ASSERT( NULL != other_s );
            result = -1;
            }
        else{
            AUTO_ASSERT( NULL == other_s );
            result = 1;
            }
        }
    }
return result;
}

uint64_t np02_boundary::hash( const uint64_t& h_in ) const{
uint64_t h = h_in;
#if defined( CF01_SUPPORT )
h = cf01_obj_hash( h, m_alloc_idx );
#else
h += reinterpret_cast<cf01_uint64>(m_alloc_idx);
h ^= ((h << 27) | (h >> 37));
#endif
if( NULL == m_segs_head ){
    AUTO_ASSERT(NULL == m_segs_tail);
    }
else{
    AUTO_ASSERT(NULL != m_segs_tail);
    h = m_segs_head->hash(h);
    const np02_boundary_seg *b = m_segs_head->get_next();
    if( NULL != b){
        const np02_boundary_seg * const b_end = (m_segs_head == m_segs_tail) ?
            m_segs_tail : NULL;
        while( b_end != b ){
            h = b->hash(h);
            b = b->get_next();
            }
        }
    }
return h;
}

double np02_boundary::get_distance_from_xy( const np02_xy& xy,
    np02_xy *near_xy, double *d_from2, bool *valid )const{
return get_distance_from_xy_exh_srch( xy, near_xy, d_from2, valid );
}

/* exhaustive search */
double np02_boundary::get_distance_from_xy_exh_srch( const np02_xy& xy,
    np02_xy *near_xy, double *d_from2, bool *valid )const{
double d_from = std::numeric_limits<double>::max();
np02_xy nrxy(0.0,0.0);
double df2 = std::numeric_limits<double>::max();
bool v = false;

np02_xy b_nrxy(0.0, 0.0);
bool b_valid = false;
double b_d = 0.0;
double b_d2 = 0.0;
bool better_answer_found = false;
if( NULL == m_segs_head ){
    AUTO_ASSERT(NULL == m_segs_tail);
    }
else{
    AUTO_ASSERT(NULL != m_segs_tail);
    b_d = m_segs_head->get_distance_from_xy( xy, &b_nrxy, &b_d2, &b_valid );
    if( b_valid ){
        v = true;
        better_answer_found = false;
        if( fabs(b_d) < fabs(d_from) ){
            better_answer_found = true;
            }
        else if( fabs(b_d) == fabs(d_from) ){
            if( b_d2 < df2 ){
                better_answer_found = true;
                }
            else if( ( b_d2 == df2 ) &&
                ( b_nrxy < nrxy ) ){
                better_answer_found = true;
                }
            }
        if( better_answer_found ){
            d_from = b_d;
            nrxy = b_nrxy;
            df2 = b_d2;
            }
        }
    const np02_boundary_seg *b = m_segs_head->get_next();
    if( NULL != b){
        const np02_boundary_seg * const b_end = (m_segs_head == m_segs_tail) ?
            m_segs_tail : NULL;
        while( b_end != b ){
            b_d = b->get_distance_from_xy( xy, &b_nrxy, &b_d2, &b_valid );
            if( b_valid ){
                v = true;
                better_answer_found = false;
                if( fabs(b_d) < fabs(d_from) ){
                    better_answer_found = true;
                    }
                else if( fabs(b_d) == fabs(d_from) ){
                    if( b_d2 < df2 ){
                        better_answer_found = true;
                        }
                    else if( ( b_d2 == df2 ) &&
                        ( b_nrxy < nrxy ) ){
                        better_answer_found = true;
                        }
                    }
                if( better_answer_found ){
                    d_from = b_d;
                    nrxy = b_nrxy;
                    df2 = b_d2;
                    }
                }
            b = b->get_next();
            }
        }
    }

if( NULL != near_xy ){
    *near_xy = nrxy;
    }
if( NULL != valid ){
    *valid = v;
    }
if( NULL != d_from2 ){
    *d_from2 = df2;
    }
return d_from;
}

int np02_boundary::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const{
int err_cnt = 0;

if( NULL != m_shp_alloc ){
    const np02_boundary * const sa_bdry = 
        m_shp_alloc->alloc_get_boundary_by_idx(m_alloc_idx);
    if( this != sa_bdry ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary: this=%llx  != (shp_alloc=%llx)->boundary_by_idx(i=%i)=%llx\n",
            this, m_shp_alloc, m_alloc_idx, sa_bdry );
        }
    }

if( NULL != m_owner ){
    const np02_shp_alloc *ownr_sa = m_owner->get_shp_alloc();
    if( ownr_sa != m_shp_alloc ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary: this=%llx  shp_alloc=%llx  !=  (owner=%llx)->shp_alloc=%llx\n",
            this, m_shp_alloc, m_owner, ownr_sa );
        }
    }

err_cnt += verify_data_segs_head_tail(err_msg, err_msg_capacity, err_msg_pos);

const np02_boundary_seg *s = m_segs_head;
while ( NULL != s ){
    err_cnt += s->verify_data(err_msg, err_msg_capacity, err_msg_pos);
    if( this != s->get_owner() ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary: this=%llx  !=  (boundary_seg=%llx)->owner=%llx\n",
            this, s, s->get_owner());
        }
    if( m_shp_alloc != s->get_shp_alloc() ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary: this=%llx shp_alloc=%llx  !=  "
            "(boundary_seg=%llx)->shp_alloc=%llx\n",
            this, m_shp_alloc, s, s->get_shp_alloc());
        }
    s = advance_boundary_seg_citr(s);
    }

if( NULL != m_rgn_bdry_prev ){
    if( this != m_rgn_bdry_prev->m_rgn_bdry_next ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary: this=%llx  !=  (rgn_bdry_prev=%llx)->rgn_bdry_next=%llx\n",
            this, m_rgn_bdry_prev, m_rgn_bdry_prev->m_rgn_bdry_next );
        }
    if( m_owner != m_rgn_bdry_prev->m_owner ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary:(this=%llx)->owner=%llx  !=  (rgn_bdry_prev=%llx)->owner=%llx\n",
            this, m_owner, m_rgn_bdry_prev, m_rgn_bdry_prev->m_owner );
        }
    if( m_shp_alloc != m_rgn_bdry_prev->m_shp_alloc ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary:(this=%llx)->shp_alloc=%llx  !=  "
            "(rgn_bdry_prev=%llx)->shp_alloc=%llx\n",
            this, m_shp_alloc, m_rgn_bdry_prev, m_rgn_bdry_prev->m_shp_alloc );
        }
    }

if( NULL != m_rgn_bdry_next ){
    if( this != m_rgn_bdry_next->m_rgn_bdry_prev ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary: this=%llx  !=  (rgn_bdry_next=%llx)->rgn_bdry_prerv=%llx\n",
            this, m_rgn_bdry_next, m_rgn_bdry_next->m_rgn_bdry_prev );
        }
    if( m_owner != m_rgn_bdry_next->m_owner ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary:(this=%llx)->owner=%llx  !=  (rgn_bdry_next=%llx)->owner=%llx\n",
            this, m_owner, m_rgn_bdry_next, m_rgn_bdry_next->m_owner );
        }
    if( m_shp_alloc != m_rgn_bdry_next->m_shp_alloc ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary:(this=%llx)->shp_alloc=%llx  !=  "
            "(rgn_bdry_next=%llx)->shp_alloc=%llx\n",
            this, m_shp_alloc, m_rgn_bdry_next, m_rgn_bdry_next->m_shp_alloc );
        }
    }

return err_cnt;
}

int np02_boundary::verify_data_segs_head_tail( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;
if( NULL == m_segs_head ){
    if ( NULL != m_segs_tail ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary: this=%llx  segs_head NULL  segs_tail=%llx\n",
            this, m_segs_tail );
        }
    } 
else if ( NULL == m_segs_tail ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "boundary: this=%llx  segs_head=%llx  segs_tail NULL\n",
         this, m_segs_head );
    }
else{
    if( m_segs_head == m_segs_tail ){
        if( m_segs_head->get_prev() == NULL){
            if( m_segs_head->get_next() == NULL){
                /* single element open (non-loop) sequence OK */
                }
            else{
                ++err_cnt;
                np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                    "boundary: this=%llx  segs_head=segs_tail=%llx  "
                    "prev NULL   next=%llx\n",
                     this, m_segs_head, m_segs_head->get_next() );
                }
            } 
        else if(m_segs_head->get_next() == NULL){
                ++err_cnt;
                np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                    "boundary: this=%llx  segs_head=segs_tail=%llx  "
                    "prev=%llx   next NULL\n",
                     this, m_segs_head, m_segs_head->get_prev() );
            }
        else{
            /* loop OK */
            }
        } 
    else{
        /* open (non-loop) sequence, size > 1 */
        if(m_segs_head->get_prev() != NULL ){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "boundary: this=%llx  segs_head=%llx != segs_tail=%llx  "
                "segs_head->prev=%llx\n",
                 this, m_segs_head, m_segs_tail, m_segs_head->get_prev() );
            }
        if(m_segs_head->get_next() == NULL ){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "boundary: this=%llx  segs_head=%llx != segs_tail=%llx  "
                "segs_head->next NULL\n",
                 this, m_segs_head, m_segs_tail );
            }
        if(m_segs_tail->get_prev() == NULL ){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "boundary: this=%llx  segs_head=%llx != segs_tail=%llx  "
                "segs_tail->prev NULL\n",
                 this, m_segs_head, m_segs_tail );
            }
        if(m_segs_tail->get_next() != NULL ){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "boundary: this=%llx  segs_head=%llx != segs_tail=%llx  "
                "segs_tail->next=%llx\n",
                 this, m_segs_head, m_segs_tail, m_segs_tail->get_next() );
            }
        }
    }
return err_cnt;
}

int np02_boundary::verify_loc_grid_properly_contains_boundary(
    const np02_loc_grid* loc_grid, const np02_boundary* boundary,
    char* err_msg, const size_t err_msg_capacity,
    size_t* err_msg_pos){
int err_cnt = 0;
if (NULL == boundary) {
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "boundary NULL\n" );
    }
else {
    const np02_boundary_seg* b = boundary->get_boundary_segs_head();
    while(NULL != b) {
        err_cnt += np02_boundary_seg::verify_loc_grid_properly_contains_boundary_seg(
            loc_grid, b, err_msg, err_msg_capacity, err_msg_pos);

        b = boundary->advance_boundary_seg_citr(b);
        }
    }
return err_cnt;
}

std::ostream& np02_boundary::ostream_output(std::ostream& os) const{
/* TODO: implement */
return os;
}

void np02_boundary::write_bmp_file(const np02_xy& xy_min,
    const double& pixel_num, const np02_bmp_color& color,
    np02_bmp_file *bmp_file) const{
/* TODO: implement */
}

void np02_boundary::write_bmp_file_edge(const np02_xy& xy_min,
    const double& pixel_num, const np02_bmp_color& color,
    np02_bmp_file *bmp_file) const{
/* TODO: implement */
}

void np02_boundary::write_dxf_file(const std::string& layer,
    const uint8_t& color, np02_dxf_file *dxf_file) const{
const np02_boundary_seg *s = m_segs_head;
while ( NULL != s ){
    s->write_dxf_file(layer, color, dxf_file);
    s = advance_boundary_seg_citr(s);
    }
}

/*
Initialize boundary as a n-sided star

                     +                          +
                    /:\                         :
                   / : \                        : \
                   / : \                        :alpha
                  /  :R \                       :   \
    +--------------------+------------+         :
       \         /   :  .\         /            :     \
          \     /    :.r  \     /               :
             \ /     +     \ /                  :       \
              / \   ctr   / \                  R:
              /    \   /    \                   :         \
             /      / \      \                  :
            /    /       \    \                 :     gamma+
            / /             \ \                 :         / 
           +                   +                :       /   
                                                :beta /r    
                                                :   /       
                                                : /         
                                                +           
                                               ctr          

     n                      n=5                           
alpha = 90/n                18 deg
beta =  180/n               36
gamma = 180-(270/n)         126
gamma = 180-(3*alpha)
r/sin(alpha) = R/sin(gamma)
r=R*(sin(alpha)/sin(gamma))
*/
void np02_boundary::init_star(const np02_xy& ctr, const double& radius,
    const size_t& n){
clear_boundary_segs();
if( (radius > 0.0) && (n > 2) ){
    /* small radius */
    const double alpha_deg = 90.0 / static_cast<double>(n);
    const double gamma_deg = 180.0 - ( 3.0 * alpha_deg );
    static const double rad_per_deg = 3.1415926535897932384626433832795/180.0;
    const double alpha_rad = rad_per_deg * alpha_deg;
    const double gamma_rad = rad_per_deg * gamma_deg;
    const double sin_alpha = sin( alpha_rad );
    const double sin_gamma = sin( gamma_rad );
    AUTO_ASSERT( alpha_deg > 0.0 );
    AUTO_ASSERT( alpha_deg < 180.0 );
    AUTO_ASSERT( sin_alpha > 0.0 );
    AUTO_ASSERT( gamma_deg > 0.0 );
    AUTO_ASSERT( gamma_deg < 180.0 );
    AUTO_ASSERT( sin_gamma > 0.0 );
    const double r_small_ratio = sin_alpha/sin_gamma;
    const double r_small = radius * r_small_ratio;

    /* construct edges */
    const double beta_rad = 2.0 * alpha_rad;
    const np02_xy p_end(ctr.get_x(), ctr.get_y() + radius);
    np02_xy p0(0.0,0.0);
    np02_xy p1 = p_end;
    const size_t i_end = 2 * n;
    size_t i = 0;
    while( i < i_end ){
        ++i;
        p0 = p1;
        if( i < i_end ){
            const double theta_rad = static_cast<double>(i) * beta_rad;
            const double cos_theta = cos(theta_rad);
            const double sin_theta = sin(theta_rad);
            const double& r = ( 0 == ( i % 2 ) ) ? radius : r_small;
            p1.set_x( ctr.get_x() - ( r * sin_theta ) );
            p1.set_y( ctr.get_y() + ( r * cos_theta ) );
            }
        else{
            p1 = p_end;
            }
        np02_line_seg *edge = (NULL == m_shp_alloc) ?
            new np02_line_seg() : m_shp_alloc->alloc_line_seg();
        edge->init(p0, p1, 0.0);
        np02_boundary_seg *b_seg = (NULL == m_shp_alloc) ?
            new np02_boundary_seg() : m_shp_alloc->alloc_boundary_seg();
        b_seg->set_shape(edge);
        add_boundary_seg_loop(b_seg);
        }
    }
}

void np02_boundary::init_rectangle(const np02_xy& rect_xy_min,
    const np02_xy& rect_xy_max){
clear_boundary_segs();
if( ( rect_xy_min.get_x() < rect_xy_max.get_x() ) &&
    ( rect_xy_min.get_y() < rect_xy_max.get_y() ) ){
    np02_xy p0(0.0,0.0);
    np02_xy p1 = rect_xy_min;
    for( size_t i = 0; i < 4; ++i ){
        p0 = p1;
        switch( i ){
        case 0:
            p1.set_x( rect_xy_max.get_x() );
            p1.set_y( rect_xy_min.get_y() );
            break;
        case 1:
            p1 = rect_xy_max;
            break;
        case 2:
            p1.set_x( rect_xy_min.get_x() );
            p1.set_y( rect_xy_max.get_y() );
            break;
        case 3:
        default:
            p1 = rect_xy_min;
            break;
        }
        np02_line_seg *edge = (NULL == m_shp_alloc) ?
            new np02_line_seg() : m_shp_alloc->alloc_line_seg();
        edge->init(p0, p1, 0.0);
        np02_boundary_seg *b_seg = (NULL == m_shp_alloc) ?
            new np02_boundary_seg() : m_shp_alloc->alloc_boundary_seg();
        b_seg->set_shape(edge);
        add_boundary_seg_loop(b_seg);
        }
    }
}


np02_region::np02_region(): m_shp_alloc(NULL), m_alloc_idx(0), m_owner(NULL),
    m_free_chain_next(NULL), m_boundaries_head(NULL), m_boundaries_tail(NULL){
}

np02_region::~np02_region(){
destruct();
}

void np02_region::destruct(){
clear_boundaries();
}

void np02_region::clear_boundaries(){
AUTO_ASSERT( 0 == verify_data_boundaries_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ))
np02_boundary *b = m_boundaries_head;
while ( NULL != b ) {
    AUTO_ASSERT( b->get_shp_alloc() == m_shp_alloc );
    AUTO_ASSERT( b->get_owner() == this );
    np02_boundary *n = b->get_rgn_bdry_next();
    AUTO_ASSERT( (NULL == n) || (n->get_rgn_bdry_prev() == b) );
    AUTO_ASSERT( (NULL == n) == (m_boundaries_tail == b) );
    b->set_owner(NULL);
    b->set_rgn_bdry_prev(NULL);
    b->set_rgn_bdry_next(NULL);
    if( NULL == m_shp_alloc ){
        delete b;
        }
    else{
        m_shp_alloc->free_boundary(b);
        }
    b = n;
    }
m_boundaries_head = NULL;
m_boundaries_tail = NULL;
AUTO_ASSERT( 0 == verify_data_boundaries_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ))
}

void np02_region::add_boundary(np02_boundary *b){
AUTO_ASSERT(0 == verify_data_boundaries_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR()));
AA_ALWAYS_ASSERT(NULL != b)
AUTO_ASSERT( b->get_shp_alloc() == m_shp_alloc );
np02_region *ownr = b->get_owner();
if( this == ownr ){
    /* already added -- do nothing*/
    AUTO_ASSERT( NULL != m_boundaries_head );
    AUTO_ASSERT( NULL != m_boundaries_tail );
    }
else {
    if(NULL != ownr){
        ownr->remove_boundary(b);
        }
    AUTO_ASSERT( NULL == b->get_rgn_bdry_prev() );
    AUTO_ASSERT( NULL == b->get_rgn_bdry_next() );
    b->set_owner(this);
    if(NULL == m_boundaries_tail){
        AUTO_ASSERT( NULL == m_boundaries_head );
        m_boundaries_head = b;
        m_boundaries_tail = b;
        }
    else{
        AUTO_ASSERT( NULL != m_boundaries_head );
        AUTO_ASSERT( NULL == m_boundaries_head->get_rgn_bdry_prev() );
        AUTO_ASSERT( NULL == m_boundaries_tail->get_rgn_bdry_next() );
        m_boundaries_tail->set_rgn_bdry_next(b);
        b->set_rgn_bdry_prev(m_boundaries_tail);
        m_boundaries_tail = b;
        }
    }
AUTO_ASSERT( 0 == verify_data_boundaries_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ));
}

void np02_region::remove_boundary(np02_boundary *b){
AUTO_ASSERT( 0 == verify_data_boundaries_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ))
AUTO_ASSERT(NULL != b);
AUTO_ASSERT((NULL == b) ||this == b->get_owner());
if( (NULL != b) && (this == b->get_owner()) ){
    b->set_owner(NULL);
    np02_boundary *p = b->get_rgn_bdry_prev();
    np02_boundary *n = b->get_rgn_bdry_next();
    if( NULL == p ){
        AUTO_ASSERT( b == m_boundaries_head );
        m_boundaries_head = n;
        } 
    else{
        AUTO_ASSERT( b != m_boundaries_head );
        p->set_rgn_bdry_next(n);
        b->set_rgn_bdry_prev(NULL);
        }
    if( NULL == n ){
        AUTO_ASSERT( b == m_boundaries_tail );
        m_boundaries_tail = p;
        } 
    else{
        AUTO_ASSERT( b != m_boundaries_tail );
        n->set_rgn_bdry_prev(p);
        b->set_rgn_bdry_next(NULL);
        }
    AUTO_ASSERT(NULL == b->get_rgn_bdry_prev());
    AUTO_ASSERT(NULL == b->get_rgn_bdry_next());
    }
AUTO_ASSERT( 0 == verify_data_boundaries_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ))
}

size_t np02_region::get_boundary_count() const{
size_t c = 0;
const np02_boundary *b = m_boundaries_head;
while( NULL != b ){
    ++c;
    b = b->get_rgn_bdry_next();
    }
return c;
}

size_t np02_region::get_shape_count() const{
size_t c = 0;
const np02_boundary *b = m_boundaries_head;
while( NULL != b ){
    c += b->get_shape_count();
    b = b->get_rgn_bdry_next();
    }
return c;
}

/* stop counting when shape count reaches threshold */
void np02_region::get_shape_count_lower_bound(
    const size_t& max_shape_count_l_b, size_t *count_l_b ) const{
if( NULL != count_l_b ){
    const np02_boundary *b = m_boundaries_head;
    while( ( NULL != b ) && ( *count_l_b < max_shape_count_l_b ) ){
        b->get_shape_count_lower_bound(max_shape_count_l_b, count_l_b);
        b = b->get_rgn_bdry_next();
        }
    }
}

np02_shape *np02_region::get_first_shape() const{
np02_shape *shape = NULL;
const np02_boundary *b = m_boundaries_head;
while( (NULL != b) && (NULL == shape) ){
    shape = b->get_first_shape();
    b = b->get_rgn_bdry_next();
    }
return shape;
}

np02_loc_grid *np02_region::get_loc_grid() const{
np02_shape *shape = get_first_shape();
np02_loc_grid *loc_grid = ( NULL == shape ) ? NULL : shape->get_loc_grid();
return loc_grid;
}

void np02_region::insert_region_in_loc_grid(np02_loc_grid *loc_grid){
np02_boundary *b = m_boundaries_head;
while( NULL != b ){
    b->insert_boundary_in_loc_grid(loc_grid);
    b = b->get_rgn_bdry_next();
    }
AUTO_ASSERT( ( NULL == get_first_shape() ) || ( get_loc_grid() == loc_grid ) );
}

void np02_region::remove_region_from_loc_grid(){
np02_boundary *b = m_boundaries_head;
while( NULL != b ){
    b->remove_boundary_from_loc_grid();
    b = b->get_rgn_bdry_next();
    }
AUTO_ASSERT( get_loc_grid() == NULL );
}

void np02_region::get_bb(np02_xy *xy_min, np02_xy *xy_max, bool *valid) const{
np02_xy xy_low(0.0,0.0);
np02_xy xy_high(0.0,0.0);
np02_xy b_xy_low(0.0,0.0);
np02_xy b_xy_high(0.0,0.0);
bool ok = false;
bool b_ok = false;
const np02_boundary *b = m_boundaries_head;
while( NULL != b ){
    if( ok ){
        b->get_bb(&b_xy_low,&b_xy_high,&b_ok);
        if( b_ok ){
            np02_shape::bb_merge(xy_low, xy_high, b_xy_low, b_xy_high,
                &xy_low, &xy_high);
            }
        }
    else{
        b->get_bb(&xy_low,&xy_high,&ok);
        }
    b = b->get_rgn_bdry_next();
    }

AUTO_ASSERT(xy_low.get_x() <= xy_high.get_x());
AUTO_ASSERT(xy_low.get_y() <= xy_high.get_y());

if( NULL != xy_min ){
    *xy_min = xy_low;
    }
if( NULL != xy_max ){
    *xy_max = xy_high;
    }
if( NULL != valid ){
    *valid = ok;
    }
}

double np02_region::get_boundary_gap_sum() const{
double gap_sum = 0.0;
const np02_boundary *b = m_boundaries_head;
while( NULL != b ){
    gap_sum += b->get_boundary_gap_sum();
    b = b->get_rgn_bdry_next();
    }
return gap_sum;
}

void np02_region::translate(const np02_xy& dxy){
np02_boundary *b = m_boundaries_head;
while( NULL != b ){
    b->translate(dxy);
    b = b->get_rgn_bdry_next();
    }
}

void np02_region::translate_no_loc_grid(const np02_xy& dxy){
np02_boundary *b = m_boundaries_head;
while( NULL != b ){
    b->translate_no_loc_grid(dxy);
    b = b->get_rgn_bdry_next();
    }
}

void np02_region::rotate(const np02_xy& rot_ctr, const double& rot_deg){
np02_boundary *b = m_boundaries_head;
while( NULL != b ){
    b->rotate(rot_ctr, rot_deg);
    b = b->get_rgn_bdry_next();
    }
}

void np02_region::rotate_no_loc_grid(const np02_xy& rot_ctr,
    const double& rot_deg){
np02_boundary *b = m_boundaries_head;
while( NULL != b ){
    b->rotate_no_loc_grid(rot_ctr, rot_deg);
    b = b->get_rgn_bdry_next();
    }
}

void np02_region::invert(){
np02_boundary *b = m_boundaries_head;
while( NULL != b ){
    b->invert();
    b = b->get_rgn_bdry_next();
    }
}

np02_region *np02_region::copy_region(np02_shp_alloc *shp_alloc) const{
np02_region *r = (NULL == shp_alloc) ?
    new np02_region() : shp_alloc->alloc_region();
r->init_from_region(this);
return r;
}

void np02_region::init_from_region( const np02_region *master ){
AA_CALL_DEPTH_BLOCK();
if( (NULL != master) && (this != master) ){
    clear_boundaries();
    np02_boundary *b = master->m_boundaries_head;
    while( NULL != b ){
        np02_boundary *b_copy = b->copy_boundary(m_shp_alloc);
        add_boundary(b_copy);
        b = b->get_rgn_bdry_next();
        }
    AUTO_ASSERT( 0 == compare_region(master) );
    }
}

int np02_region::compare_region( const np02_region *other ) const{
int result = 0;
if( NULL == other ){
    result = 1;
    }
else{
    const np02_boundary *b = m_boundaries_head;
    const np02_boundary *other_b = other->m_boundaries_head;
    while( ( NULL != b ) && ( NULL != other_b ) && ( 0 == result ) ){
        result = b->compare_boundary(other_b);
        b = b->get_rgn_bdry_next();
        other_b = other_b->get_rgn_bdry_next();
        }
    if( ( 0 == result ) && ( ( ( NULL != b ) || ( NULL != other_b ) ) ) ){
        if( NULL == b ){
            AUTO_ASSERT( NULL != other_b );
            result = -1;
            }
        else{
            AUTO_ASSERT( NULL == other_b );
            result = 1;
            }
        }
    }
return result;
}

uint64_t np02_region::hash( const uint64_t& h_in ) const{
uint64_t h = h_in;
#if defined( CF01_SUPPORT )
h = cf01_obj_hash( h, m_alloc_idx );
#else
h += reinterpret_cast<cf01_uint64>(m_alloc_idx);
h ^= ((h << 27) | (h >> 37));
#endif
const np02_boundary *b = m_boundaries_head;
while ( NULL != b ) {
    h = b->hash(h);
    b = b->get_rgn_bdry_next();
    }
return h;
}

double np02_region::get_distance_from_xy(const np02_xy& xy,
   np02_xy *near_xy, double *d_from2, bool *valid )const{
np02_shape_vec temp_shape_vec;
const double d_from = get_distance_from_xy( xy, near_xy, d_from2, valid,
    &temp_shape_vec );
return d_from;
}

double np02_region::get_distance_from_xy(
    const np02_xy& xy, np02_xy *near_xy, double *d_from2, bool *valid,
    np02_shape_vec *temp_shape_vec )const{
static const double d_from_default = std::numeric_limits<double>::max();
double d_from = d_from_default;
np02_xy nr_xy(0.0,0.0);
double nr_d2=d_from_default;
bool should_do_exhaustive_search = false;
bool search_valid = false;
const np02_loc_grid *loc_grid = get_loc_grid();
if( NULL == loc_grid ){
    should_do_exhaustive_search = true;
    }
else{
    static const size_t shape_count_threshold = 32;
    size_t count_l_b = 0;
    get_shape_count_lower_bound( shape_count_threshold, &count_l_b );
    if( 0 == count_l_b ){
        AUTO_ASSERT( !should_do_exhaustive_search );
        AUTO_ASSERT( !search_valid );
        }
    else if( count_l_b < shape_count_threshold ){
        should_do_exhaustive_search = true;
        }
    else{
        /* search locator grid*/
        AA_ALWAYS_ASSERT( NULL != temp_shape_vec );

        const np02_loc_grid_dim& loc_grid_dim = loc_grid->get_loc_grid_dim();
        np02_xy search_xy_min;
        np02_xy search_xy_max;
        np02_xy near_shp_xy;
        double near_shp_d = 0.0;
        double near_shp_d2 = 0.0;
        bool better_answer_found = false;
        const double start_search_radius = loc_grid_dim.get_sq_size();
        if( start_search_radius <= 0.0 ){
            should_do_exhaustive_search = true;
            }
        const uint16_t fourth_w_h = 
            (loc_grid_dim.get_w() + loc_grid_dim.get_h())/4;
        const double end_search_radius = (fourth_w_h > 1 ) ? 
            (start_search_radius * static_cast<double>(fourth_w_h)) :
            (start_search_radius * 2.0);
        double search_radius = start_search_radius;
        while( ( !search_valid ) && ( !should_do_exhaustive_search ) &&
            ( search_radius <= end_search_radius ) ){
            /* get shapes within search radius */
            search_xy_min.set_x( xy.get_x() - search_radius );
            search_xy_min.set_y( xy.get_y() - search_radius );
            search_xy_max.set_x( xy.get_x() + search_radius );
            search_xy_max.set_y( xy.get_y() + search_radius );
            temp_shape_vec->resize(0);
            loc_grid->get_shapes_in_bb(search_xy_min, search_xy_max,
                temp_shape_vec );

            /* check distance to each shape belonging to region */
            np02_shape_vec_citr local_shape_itr = temp_shape_vec->begin();
            for(; temp_shape_vec->end() != local_shape_itr; ++local_shape_itr){
                const np02_shape *local_shape = *local_shape_itr;
                AUTO_ASSERT(NULL != local_shape);
                if( this == local_shape->get_region() ){
                    near_shp_d = local_shape->get_bseg_dist_from_xy(xy,
                        &near_shp_xy, &near_shp_d2);
                    better_answer_found = false;
                    if( fabs(near_shp_d) < fabs(d_from) ){
                        better_answer_found = true;
                        }
                    else if( fabs(near_shp_d) == fabs(d_from) ){
                        if( near_shp_d2 < nr_d2 ){
                            better_answer_found = true;
                            }
                        else if( ( near_shp_d2 == nr_d2 ) &&
                            ( near_shp_xy < nr_xy ) ){
                            better_answer_found = true;  
                            }
                        }
                    if( better_answer_found ){
                        d_from = near_shp_d;
                        nr_xy = near_shp_xy;
                        nr_d2 = near_shp_d2;
                        }
                    }
                }

            AUTO_ASSERT(!search_valid);
            AUTO_ASSERT(!should_do_exhaustive_search);
            if( fabs(d_from) < search_radius ){
                /* success */
                search_valid = true;
                }
            else if( search_radius >= end_search_radius ){
                /* failed to find near point within largest
                search radius.  Do exhaustive search */
                search_radius *= 2.0; /* exit locator grid search */
                should_do_exhaustive_search = true;
                }
            else if( d_from_default == d_from) {
                /* nothing found.  double search radius */
                search_radius *= 2.0;
                if( search_radius > end_search_radius ){
                    /* do one final locator grid search */
                    search_radius = end_search_radius;
                    }
                }
            else{
                /* expand search radius to include any possible
                closer shapes */
                search_radius = fabs(d_from) + ( start_search_radius / 4.0 );
                if (search_radius > end_search_radius) {
                    /* do one final locator grid search */
                    search_radius = end_search_radius;
                    }
                }
            }
        temp_shape_vec->resize(0);
        }
    }

if( should_do_exhaustive_search ){
    d_from = get_distance_from_xy_exh_srch( xy, &nr_xy, &nr_d2, &search_valid );
    }
else{
    AUTO_ASSERT( 0 == verify_distance_from_xy_result( xy, d_from, nr_xy, nr_d2,
        search_valid,AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),AA_ERR_BUF_POS_PTR()));
    }

if( NULL != valid ){ *valid = search_valid; }
if( NULL != near_xy ){ *near_xy = nr_xy; }
if( NULL != d_from2){ *d_from2 = nr_d2; }

return d_from;
}

/* exhaustive search */
double np02_region::get_distance_from_xy_exh_srch( const np02_xy& xy,
    np02_xy *near_xy, double *d_from2, bool *valid )const{
static const double d_from_default = std::numeric_limits<double>::max();
double d_from = d_from_default;
np02_xy nrxy(0.0,0.0);
double df2 = d_from_default;
bool better_answer_found = false;
bool v = false;
const np02_boundary *b = m_boundaries_head;
while ( NULL != b ) {
    np02_xy b_nrxy(0.0, 0.0);
    bool b_valid = false;
    double b_df2 = 0.0;
    const double d = 
        b->get_distance_from_xy_exh_srch( xy, &b_nrxy, &b_df2, &b_valid );
    if( b_valid ){
        v = true;
        better_answer_found = false;
        if (fabs(d) < fabs(d_from)) {
            better_answer_found = true;
            }
        else if(fabs(d) == fabs(d_from)) {
            if (b_df2 < df2) {
                better_answer_found = true;
                }
            else if ((b_df2 == df2) && (b_nrxy < nrxy)){
                better_answer_found = true;
                }
            }
        if (better_answer_found) {
            d_from = d;
            nrxy = b_nrxy;
            df2 = b_df2;
        }
    }
    b = b->get_rgn_bdry_next();
    }
if( NULL != near_xy ){
    *near_xy = nrxy;
    }
if (NULL != valid) {
    *valid = v;
}
if (NULL != d_from2) {
    *d_from2 = df2;
}
return d_from;
}

/* result = 0 => */
int np02_region::verify_distance_from_xy_result( const np02_xy& xy, 
    const double& d_from, np02_xy near_xy, const double& d_from2,
    const bool& valid, char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const{
int err_cnt = 0;
np02_xy nr_xy(0.0,0.0);
bool vld = false;
double df2 = 0.0;
const double df = get_distance_from_xy_exh_srch(xy, &nr_xy, &df2, &vld);
if((df != d_from) || (nr_xy != near_xy) || (df2 != d_from2) || (vld != valid)){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "region: this=%llx  distance from xy (%g,%g) result error  "
        "d_from=%g<-->%g (error=%g)   near_xy=(%g,%g)<-->(%g,%g) (error=%g,%g)"
        "d_from2=%g<-->%g (error=%g)   valid=%i<-->%i\n",
        this, xy.get_x(), xy.get_y(), d_from, df, d_from - df,
        near_xy.get_x(), near_xy.get_y(), nr_xy.get_x(), nr_xy.get_y(), 
        near_xy.get_x() - nr_xy.get_x(), near_xy.get_y() - nr_xy.get_y(),
        d_from2, df2, d_from2 - df2,
        static_cast<int>(valid), static_cast<int>(vld) );
    }
return err_cnt;
}

int np02_region::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const{
int err_cnt = 0;

if( NULL != m_shp_alloc ){
    const np02_region * sa_rgn =
        m_shp_alloc->alloc_get_region_by_idx(m_alloc_idx);
    if( this != sa_rgn ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "region: this=%llx  !=  "
            "(shp_alloc=%llx)->get_region_by_idx(i=%i)=%llx\n",
            this, m_shp_alloc, m_alloc_idx, sa_rgn );
        }
    }

if( NULL != m_owner ){
    const np02_region * ownr_rgn = m_owner->get_region();
    if( this != ownr_rgn ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "region: this=%llx  !=  (owner=%llx)->region=%llx\n",
            this, m_owner, ownr_rgn );
        }
    }

const np02_boundary *b = m_boundaries_head;
while ( NULL != b ) {
    err_cnt += b->verify_data( err_msg, err_msg_capacity, err_msg_pos );
    if( this != b->get_owner() ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "region: this=%llx  !=  (boundary=%llx)->owner=%llx\n",
            this, b, b->get_owner() );
        }
    if( m_shp_alloc != b->get_shp_alloc() ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "region: (this=%llx) shp_alloc=%llx !=  (boundary=%llx)->owner=%llx\n",
            this, m_shp_alloc, b, b->get_shp_alloc() );
        }
    b = b->get_rgn_bdry_next();
    }

err_cnt += verify_data_boundaries_head_tail( err_msg, err_msg_capacity, err_msg_pos );

return err_cnt;
}


int np02_region::verify_data_boundaries_head_tail( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt = 0;
if( NULL == m_boundaries_head ){
    if ( NULL != m_boundaries_tail ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "region: this=%llx  boundaries_head NULL  boundaries_tail=%llx\n",
            this, m_boundaries_tail );
        }
    } 
else if ( NULL == m_boundaries_tail ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "region: this=%llx  boundaries_head=%llx  boundaries_tail NULL\n",
         this, m_boundaries_head );
    }
else{
    /* open (non-loop) sequence */
    if(m_boundaries_head->get_rgn_bdry_prev() != NULL ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "region: this=%llx  boundaries_head=%llx  boundaries_tail=%llx  "
            "boundaries_head->rgn_bdry_prev = %llx != NULL\n",
             this, m_boundaries_head, m_boundaries_tail,
             m_boundaries_head->get_rgn_bdry_prev() );
        }
    if (m_boundaries_tail->get_rgn_bdry_next() != NULL) {
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "region: this=%llx  boundaries_head=%llx boundaries_tail=%llx  "
            "boundaries_tail->rgn_bdry_next=%llx != NULL\n",
            this, m_boundaries_head, m_boundaries_tail,
            m_boundaries_tail->get_rgn_bdry_next());
    }
    if(m_boundaries_head != m_boundaries_tail){
        /* open (non-loop) sequence, size > 1 */
        if(m_boundaries_head->get_rgn_bdry_next() == NULL ){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "region: this=%llx  boundaries_head=%llx  boundaries_tail=%llx  "
                "boundaries_head->rgn_bdry_next NULL\n",
                 this, m_boundaries_head, m_boundaries_tail );
            }
        if(m_boundaries_tail->get_rgn_bdry_prev() == NULL ){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "region: this=%llx  boundaries_head=%llx != boundaries_tail=%llx  "
                "boundaries_tail->rgn_bdry_prev NULL\n",
                 this, m_boundaries_head, m_boundaries_tail );
            }
        }
    }

return err_cnt;
}

int np02_region::verify_loc_grid_properly_contains_region(
    const np02_loc_grid* loc_grid, const np02_region* region,
    char* err_msg, const size_t err_msg_capacity,
    size_t* err_msg_pos){
int err_cnt = 0;
if (NULL == region) {
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "region NULL\n");
}
else {
    const np02_boundary* b = region->get_boundaries_head();
    while (NULL != b) {
        err_cnt += np02_boundary::verify_loc_grid_properly_contains_boundary(
            loc_grid, b, err_msg, err_msg_capacity, err_msg_pos);
        b = b->get_rgn_bdry_next();
        }
    }
return err_cnt;
}

std::ostream& np02_region::ostream_output(std::ostream& os) const{
/* TODO: implement */
return os;
}

void np02_region::write_bmp_file(const np02_xy& xy_min,
    const double& pixel_num, const np02_bmp_color& color,
    np02_bmp_file *bmp_file, const time_t time_limit) const{
static const bool edges_only = false;
write_bmp_file_impl( xy_min, pixel_num, color, edges_only, bmp_file,
    time_limit);
}

void np02_region::write_bmp_file_edge(const np02_xy& xy_min,
    const double& pixel_num, const np02_bmp_color& color,
    np02_bmp_file *bmp_file, const time_t time_limit) const{
static const bool edges_only = true;
write_bmp_file_impl( xy_min, pixel_num, color, edges_only, bmp_file,
    time_limit);
}

void np02_region::write_dxf_file(const std::string& layer,
        const uint8_t& color, np02_dxf_file *dxf_file) const{
const np02_boundary *b = m_boundaries_head;
while ( NULL != b ) {
    b->write_dxf_file(layer, color, dxf_file);
    b = b->get_rgn_bdry_next();
    }
}

/*                canton
        +-----------fly---------+

    +   +-----------------------+-------------------------------------+     +
    |   | *   *   *   *   *   * |__________________red________________|     |
    |   |   *   *   *   *   *   |_________________white_______________|     |
    |   | *   *   *   *   *   * |                  red                |     |
 canton |   *   *   *   *   *   |-------------------------------------+     |
  hoist | *   *   *   *   *   * |                 white               |     |
    |   |   *   *   *   *   *   |-------------------------------------+     |
    |   | *   *   *   *   *   * |__________________red________________|     |
    |   |   *   *   *   *   *   |_________________white_______________|     |
    |   | *   *   *   *   *   * |                  red                |   hoist
    +   +-----------------------+-------------------------------------+     |
        |_________________________________________white_______________|     |
        |                                          red                |     |
        +-------------------------------------------------------------+     |
        |_________________________________________white_______________|     |
        |                                          red                |     |
        +-------------------------------------------------------------+     |
        |_________________________________________white_______________|     |
        |                                          red                |     |
        +-------------------------------------------------------------+     +

        +----------------------------fly------------------------------+

https://en.wikisource.org/wiki/Executive_Order_10834
*/
void np02_region::init_usa_flag( const np02_xy& bb_xy_min,
    const np02_xy& bb_xy_max, np02_region *red_rgn,
    np02_region *white_rgn, np02_region *blue_rgn ){
AUTO_ASSERT(bb_xy_min.get_x() < bb_xy_max.get_x());
AUTO_ASSERT(bb_xy_min.get_y() < bb_xy_max.get_y());
AUTO_ASSERT(NULL != red_rgn);
AUTO_ASSERT(NULL != white_rgn);
AUTO_ASSERT(NULL != blue_rgn);

static const double fly_hoist_ratio = 1.9;
static const double canton_fly_ratio = 0.4;
static const double stripe_count = 13.0;
static const double canton_hoist_ratio = 7.0 / stripe_count;
static const double star_diameter_ratio = 0.0616;
const double star_x_margin_ratio = 1.0 / 12.0;
const double star_x_spacing_ratio = 1.0 / 12.0;
const double star_y_margin_ratio = 1.0 / 12.0;
const double star_y_spacing_ratio = 1.0 / 12.0;

const double fly_a = bb_xy_max.get_x() - bb_xy_min.get_x();
const double hoist_b = bb_xy_max.get_y() - bb_xy_min.get_y();
const double fly_b = hoist_b * fly_hoist_ratio;
double hoist = 0.0;
double fly = 0.0;
if( fly_a < fly_b ){
    fly = fly_a;
    hoist = fly_a / fly_hoist_ratio;
    }
else{
    fly = fly_b;
    hoist = hoist_b;
    }
const double canton_fly = fly * canton_fly_ratio;
const double canton_hoist = hoist * canton_hoist_ratio;
const double stripe_width = hoist / stripe_count;
const double star_diameter = hoist * star_diameter_ratio;
const double star_radius = star_diameter / 2.0;
const double star_x_margin = canton_fly * star_x_margin_ratio;
const double star_x_spacing = canton_fly * star_x_spacing_ratio;
const double star_y_margin = canton_hoist * star_y_margin_ratio;
const double star_y_spacing = canton_hoist * star_y_spacing_ratio;

const np02_xy ctr = bb_xy_min.get_midpoint_to(bb_xy_max);
const np02_xy flag_xy_min(ctr.get_x()-(fly/2.0), ctr.get_y()-(hoist/2.0));
const np02_xy flag_xy_max(ctr.get_x()+(fly/2.0), ctr.get_y()+(hoist/2.0));
const np02_xy canton_xy_min( flag_xy_min.get_x(),
    flag_xy_min.get_y() + (6.0 * stripe_width) );
const np02_xy canton_xy_max( flag_xy_min.get_x() + canton_fly,
    ctr.get_y()+(hoist/2.0 ) );

if( (NULL != red_rgn) && (NULL != white_rgn) && (NULL != blue_rgn) ){
    red_rgn->clear_boundaries();
    white_rgn->clear_boundaries();
    blue_rgn->clear_boundaries();

    /* canton */
    np02_shp_alloc *blue_shp_alloc = blue_rgn->get_shp_alloc();
    np02_boundary *canton_bdry = ( NULL == blue_shp_alloc ) ?
        new np02_boundary() : blue_shp_alloc->alloc_boundary();
    canton_bdry->init_rectangle( canton_xy_min, canton_xy_max );
    blue_rgn->add_boundary(canton_bdry);

    /* stars */
    np02_shp_alloc *white_shp_alloc = white_rgn->get_shp_alloc();
    np02_boundary *white_star_0 = ( NULL == white_shp_alloc ) ?
        new np02_boundary() : white_shp_alloc->alloc_boundary();
    const np02_xy star_0_ctr( canton_xy_min.get_x() + star_x_margin,
        canton_xy_min.get_y() + star_y_margin );
    white_star_0->init_star( star_0_ctr, star_radius, 5 );
    white_rgn->add_boundary(white_star_0);
    np02_boundary *blue_hole_0 = white_star_0->copy_boundary(blue_shp_alloc);
    blue_hole_0->invert();
    blue_rgn->add_boundary(blue_hole_0);
    np02_xy dxy(0.0,0.0);
    for( size_t i = 0; i < 11; ++i ){
        dxy.set_x( static_cast<double>(i) * star_x_spacing);
        for( size_t j = 0; j < 11; ++j ){
            if ( ((i>0)||(j>0)) && (((i+j)%2)==0) ) {
                dxy.set_y( static_cast<double>(j) * star_y_spacing);
                np02_boundary *white_star =
                    white_star_0->copy_boundary(white_shp_alloc);
                white_star->translate(dxy);
                white_rgn->add_boundary(white_star);
                np02_boundary *blue_hole =
                    blue_hole_0->copy_boundary(blue_shp_alloc);
                blue_hole->translate(dxy);
                blue_rgn->add_boundary(blue_hole);
                }
            }
        }

    /* stripes */
    double stripe_y_max = flag_xy_min.get_y();
    np02_xy stripe_xy_min = flag_xy_min;
    np02_xy stripe_xy_max = flag_xy_max;
    for( size_t stripe_idx = 0; stripe_idx < 13; ++stripe_idx ){
        stripe_xy_min.set_y(stripe_y_max);
        switch( stripe_idx ){
        case 5:
            stripe_y_max = canton_xy_min.get_y();
            break;
        case 6:
            stripe_xy_min.set_x(canton_xy_max.get_x());
            stripe_y_max += stripe_width;
            break;
        case 12:
            stripe_y_max = flag_xy_max.get_y();
            break;
        default:
            stripe_y_max += stripe_width;
            break;
            }
        stripe_xy_max.set_y(stripe_y_max);
        np02_region *region = ( 0 == ( stripe_idx % 2 ) ) ? red_rgn : white_rgn;
        np02_shp_alloc *shp_alloc = region->get_shp_alloc();
        np02_boundary *stripe = ( NULL == shp_alloc ) ?
            new np02_boundary() : shp_alloc->alloc_boundary();
        stripe->init_rectangle( stripe_xy_min, stripe_xy_max );
        region->add_boundary(stripe);
        }
    }
}

void np02_region::write_bmp_file_impl(const np02_xy& xy_min,
    const double& pixel_num, const np02_bmp_color& color,
    const bool& edges_only, np02_bmp_file *bmp_file,
    const time_t time_limit) const{

enum process_mode{
    PROCESS_MODE_SPLIT,
    PROCESS_MODE_COLOR_BOX,
    PROCESS_MODE_AVOID_COLORING_BOX,
    PROCESS_MODE_COUNT
};

time_t now = time(NULL);
time_t deadline = now + time_limit + 1;
if (deadline < now) {
    deadline = std::numeric_limits<time_t>::max();
    }

static const double outside_box_threshold_ratio = 65537.0 / 65536.0;
const double outside_box_threshold_factor= ( pixel_num > 0.0 ) ? 
    ( ( outside_box_threshold_ratio / 2.0 ) / pixel_num ) : 1.0;
const double edge_d_threshold =
    ( pixel_num > 0.0 ) ? ( 0.75 / pixel_num ) : 1.0;

/* get region bounding box */
np02_xy bb_xy_min(0.0, 0.0); 
np02_xy bb_xy_max(0.0, 0.0);
bool bb_valid = false;
get_bb( &bb_xy_min, &bb_xy_max, &bb_valid );
if( !bb_valid ){
    bb_xy_min = xy_min;
    bb_xy_max = xy_min;
    }
AUTO_ASSERT(bb_xy_min.get_x() <= bb_xy_max.get_x());
AUTO_ASSERT(bb_xy_min.get_y() <= bb_xy_max.get_y());

/* get bounding box indices in bmp coordinates (include at least 
one extra pixel on all sides ) */
const int32_t i_min_a = static_cast<int32_t>(
    (bb_xy_min.get_x() - xy_min.get_x()) * pixel_num);
const int32_t i_min = bmp_file->get_clamped_i(i_min_a - 1);

const int32_t i_max_a = static_cast<int32_t>(
    (bb_xy_max.get_x() - xy_min.get_x()) * pixel_num);
const int32_t i_max = bmp_file->get_clamped_i(i_max_a + 1);

const int32_t j_min_a = static_cast<int32_t>(
    (bb_xy_min.get_y() - xy_min.get_y()) * pixel_num);
const int32_t j_min = bmp_file->get_clamped_j(j_min_a - 1);

const int32_t j_max_a = static_cast<int32_t>(
    (bb_xy_max.get_y() - xy_min.get_y()) * pixel_num);
const int32_t j_max = bmp_file->get_clamped_j(j_max_a + 1);

AUTO_ASSERT(i_min_a <= i_max_a);
AUTO_ASSERT(j_min_a <= j_max_a);
AUTO_ASSERT(i_min <= i_max);
AUTO_ASSERT(j_min <= j_max);

typedef std::pair<int32_t, int32_t> i32_box_corner;
typedef std::pair<i32_box_corner, i32_box_corner> i32_box;
typedef std::vector<i32_box> i32_box_stack;

i32_box_stack box_stack;
box_stack.reserve(64);
i32_box box(i32_box_corner(i_min,j_min),
            i32_box_corner(i_max,j_max));
i32_box box_a;
i32_box box_b;
if( bb_valid && ( NULL != bmp_file ) ){
    box_stack.push_back(box);
    }
while( !box_stack.empty() &&
    ( ( 0 == time_limit ) || ( ( now = time(NULL) ) < deadline ) ) ){
    box = box_stack.back();
    box_stack.pop_back();

    process_mode mode = PROCESS_MODE_AVOID_COLORING_BOX;

    const int32_t& ii_min = box.first.first;
    const int32_t& jj_min = box.first.second;
    const int32_t& ii_max = box.second.first;
    const int32_t& jj_max = box.second.second;
    AUTO_ASSERT(ii_min <= ii_max);
    AUTO_ASSERT(jj_min <= jj_max);

    /* if box is too big, split in half */
    static const int32_t max_working_box_size = 32;
    const int32_t ii_range = ii_max - ii_min;
    const int32_t jj_range = jj_max - jj_min;
    if ( ( ii_range > max_working_box_size ) ||
         ( jj_range > max_working_box_size ) ){
         mode = PROCESS_MODE_SPLIT;
         }
    else{ 
        /*              outside_box_threshold_dx_near
                        +-----------+
                                      
            +-----------+-----------+
            |           |           |
            |           |           |
            |           |           |
            |           |           |
            |           |(xc,yc)    |
            +-----------+-----------+       +
            |           |           |       |
            |           |           |       |
            |           |           |       |outside_box_threshold_dy_near
            |           |           |       |
            |           |           |       |
            +-----------+-----------+       +
        */         

        /* box center */
        const double ii_ctr = static_cast<double>(ii_min + ii_max)/2.0;
        const double jj_ctr = static_cast<double>(jj_min + jj_max)/2.0;
        const double xc = xy_min.get_x() + ( ii_ctr / pixel_num );
        const double yc = xy_min.get_y() + ( jj_ctr / pixel_num );
        const np02_xy box_ctr_xy(xc, yc);

        /* search for nearest point on region boundary */
        bool near_search_valid = false;
        np02_xy near_xy(0.0,0.0);
        const double d_near = get_distance_from_xy( box_ctr_xy, &near_xy, NULL,
            &near_search_valid );

        if(!near_search_valid){
            AUTO_ASSERT(false); /* search should be valid if bb is valid */
            mode = PROCESS_MODE_AVOID_COLORING_BOX;
            }
        else{
            if( ( 0 == ii_range ) && ( 0 == jj_range ) ){
                /* single pixel */
                if( edges_only ){
                    mode = ( fabs( d_near ) < edge_d_threshold ) ?
                        PROCESS_MODE_COLOR_BOX :
                        PROCESS_MODE_AVOID_COLORING_BOX;
                    }
                else{
                    mode = ( d_near <= 0 ) ? PROCESS_MODE_COLOR_BOX :
                        PROCESS_MODE_AVOID_COLORING_BOX;
                    }
                }
            else{
                const double outside_box_threshold_dx_near = 
                    static_cast<double>(ii_range)*outside_box_threshold_factor;
                const double outside_box_threshold_dy_near = 
                    static_cast<double>(jj_range)*outside_box_threshold_factor;

                const double dx_near = near_xy.get_x() - xc;
                const double dy_near = near_xy.get_y() - yc;

                if( ( fabs(dx_near) > outside_box_threshold_dx_near ) &&
                    ( fabs(dy_near) > outside_box_threshold_dy_near ) ){
                    /* near point is clearly outside box */
                    if( edges_only ){
                        const double outside_box_edge_threshold_dx_near =
                            outside_box_threshold_dx_near + edge_d_threshold;
                        const double outside_box_edge_threshold_dy_near =
                            outside_box_threshold_dy_near + edge_d_threshold;
                        mode = 
                          ((fabs(dx_near)>outside_box_edge_threshold_dx_near)&&
                           (fabs(dy_near)>outside_box_edge_threshold_dy_near))?
                            PROCESS_MODE_AVOID_COLORING_BOX:PROCESS_MODE_SPLIT;
                        }
                    else{
                        mode = ( d_near <= 0 ) ? PROCESS_MODE_COLOR_BOX :
                            PROCESS_MODE_AVOID_COLORING_BOX;
                        }
                    }
                else{
                    /* near point is near to box or inside box */
                    mode = PROCESS_MODE_SPLIT;
                    }
                }
            }
        }

    switch( mode ) {
    case PROCESS_MODE_SPLIT:
        if( ii_range > jj_range ){
            box_a.first.first   = ii_min;
            box_a.first.second  = jj_min;
            box_a.second.first  = (ii_min + ii_max)/2;
            box_a.second.second = jj_max;
            box_stack.push_back(box_a);
            AUTO_ASSERT(box_a.first.first <= box_a.second.first);
            AUTO_ASSERT(box_a.first.second <= box_a.second.second);
    
            box_b.first.first   = box_a.second.first + 1;
            box_b.first.second  = jj_min;
            box_b.second.first  = ii_max;
            box_b.second.second = jj_max;
            box_stack.push_back(box_b);
            AUTO_ASSERT(box_b.first.first <= box_b.second.first);
            AUTO_ASSERT(box_b.first.second <= box_b.second.second);
            }
        else{
            box_a.first.first   = ii_min;
            box_a.first.second  = jj_min;
            box_a.second.first  = ii_max;
            box_a.second.second = (jj_min + jj_max)/2;
            box_stack.push_back(box_a);
            AUTO_ASSERT(box_a.first.first <= box_a.second.first);
            AUTO_ASSERT(box_a.first.second <= box_a.second.second);
    
            box_b.first.first   = ii_min;
            box_b.first.second  = box_a.second.second + 1;
            box_b.second.first  = ii_max;
            box_b.second.second = jj_max;
            box_stack.push_back(box_b);
            AUTO_ASSERT(box_b.first.first <= box_b.second.first);
            AUTO_ASSERT(box_b.first.second <= box_b.second.second);
            }
        break;
    case PROCESS_MODE_COLOR_BOX:
        bmp_file->draw_box( ii_min, jj_min, ii_max, jj_max, color);
        break;
    case PROCESS_MODE_AVOID_COLORING_BOX:
    case PROCESS_MODE_COUNT:
    default:
        break;
        }
    }
}


uint64_t np02_region_vec_hash( const np02_region_vec& v, 
    const uint64_t& h_in ){
uint64_t h = h_in;
for( np02_region_vec_citr v_itr = v.begin(); v_itr != v.end(); ++v_itr ){
    const np02_region *r = *v_itr;
    if ( NULL == r) {
        static const uint8_t zero = 0;
        h = cf01_obj_hash( h, zero );
        }
    else{
        h = r->hash(h);
        }
    }
return h;
}

void np02_region_vec_copy_regions( const np02_region_vec& v, 
    np02_region_vec *destination, np02_shp_alloc *shp_alloc){
AA_CALL_DEPTH_BLOCK();
if( NULL != destination ){
    destination->reserve( destination->size() + v.size() );
    np02_region_vec_citr region_itr = v.begin();
    for(; region_itr != v.end(); ++region_itr){
        const np02_region *master = *region_itr;
        np02_region *r = (NULL == master) ? NULL : master->copy_region(shp_alloc);
        destination->push_back(r);
        }
    AUTO_ASSERT( (destination->size() > v.size()) ||
        ( 0 == np02_region_vec_compare( v, *destination ) ) );
    }
}

int np02_region_vec_compare( const np02_region_vec& x,
    const np02_region_vec& y ){
int result = 0;
np02_region_vec_citr x_itr = x.begin();
np02_region_vec_citr y_itr = y.begin();
for(; (0==result) && (x_itr!=x.end()) && (y_itr!=y.end()); ++x_itr, ++y_itr){
    const np02_region *region_x = *x_itr;
    const np02_region *region_y = *y_itr;
    if( NULL == region_x ){
        if( NULL == region_y ){
            AUTO_ASSERT( 0 == result );
            }
        else{
            result = -1;
            }
        }
    else if( NULL == region_y ){
        result = 1;
        }
    else{
        result = region_x->compare_region(region_y);
        }
    }
if( 0 == result ){
    if( x.size() < y.size() ){ result = -1; }
    else if( y.size() < x.size() ){ result = 1; }
    }
return result;
}

void np02_region_vec_get_bb(const np02_region_vec& v, np02_xy* xy_min,
    np02_xy* xy_max, bool* valid) {
bool vld = false;
if( (NULL != xy_min) && (NULL != xy_max) ){
    np02_xy r_bb_xy_min(0.0,0.0);
    np02_xy r_bb_xy_max(0.0,0.0);
    for (np02_region_vec_citr v_itr = v.begin(); v_itr != v.end(); ++v_itr) {
        const np02_region* r = *v_itr;
        if (NULL != r) {
            bool vv = false;
            r->get_bb(&r_bb_xy_min, &r_bb_xy_max, &vv);
            if (vv) {
                if (!vld) {
                    vld = true;
                    *xy_min = r_bb_xy_min;
                    *xy_max = r_bb_xy_max;
                    }
                else {
                    np02_shape::bb_merge(*xy_min, *xy_max, r_bb_xy_min,
                        r_bb_xy_max, xy_min, xy_max);
                     }
                }
            }
        }
    }
if (NULL != valid) {
    *valid = vld;
    }
}

size_t np02_region_vec_get_shape_count(const np02_region_vec& v) {
size_t count = 0;
for (np02_region_vec_citr v_itr = v.begin(); v_itr != v.end(); ++v_itr) {
    const np02_region* r = *v_itr;
    if (NULL != r) {
        count += r->get_shape_count();
        }
    }
return count;
}

} /* namespace np02 */
