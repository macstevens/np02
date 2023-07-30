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
            "boundary_seg: this=%x  != "
            "(shp_alloc=%x) boundary_seg_by_idx(idx=%i)=%x\n",
            this, m_shp_alloc, m_alloc_idx, b_seg ); 
        }
    }

if( NULL != m_shape ){
    err_cnt += m_shape->verify_data(err_msg,err_msg_capacity,err_msg_pos);
    const np02_shape_owner *shp_ownr = m_shape->get_shape_owner();
    if( this != shp_ownr ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary_seg: this=%x  != (shape=%x) owner=%x\n",
            this, m_shape, shp_ownr ); 
        }
    const np02_shp_alloc *sa = m_shape->get_shp_alloc();
    if( m_shp_alloc != sa ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary_seg: (this=%x) shp_alloc=%x != "
            "(shape=%x) shp_alloc=%x\n",
            this, m_shp_alloc, m_shape, sa );
        }
    }

if( NULL != m_prev ){
    if( this != m_prev->m_next ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary_seg: this=%x != (prev=%x)->(next=%x)\n",
            this, m_prev, m_prev->m_next ); 
        }
    if( m_owner != m_prev->m_owner ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary_seg: (this=%x) (owner=%x) != (prev=%x)->(owner=%x)\n",
            this, m_owner, m_prev, m_prev->m_owner ); 
        }
    }

if( NULL != m_next ){
    if( this != m_next->m_prev ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary_seg: this=%x != (next=%x)->(prev=%x)\n",
            this, m_next, m_next->m_prev ); 
        }
    if( m_owner != m_next->m_owner ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary_seg: (this=%x) (owner=%x) != (next=%x)->(owner=%x)\n",
            this, m_owner, m_next, m_next->m_owner ); 
        }
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

void np02_boundary_seg::write_dxf_file(const std::string& layer,
    const uint8_t& color, np02_dxf_file *dxf_file) const{
if( NULL != m_shape ){
    m_shape->write_dxf_file(layer, color, dxf_file);
    }
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
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ))
AUTO_ASSERT( 0 == verify_data_segs_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ))
}

/* add boundary_seg, make open (non-loop) sequence */
void np02_boundary::add_boundary_seg_open(np02_boundary_seg *s){
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
AUTO_ASSERT( 0 == verify_data_segs_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ))
}

void np02_boundary::remove_boundary_seg(np02_boundary_seg *s){
AUTO_ASSERT( 0 == verify_data_segs_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ))
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
        if( NULL != n ){ n->set_next(p); }
        }
    s->set_owner(NULL);
    s->set_prev(NULL);
    s->set_next(NULL);
    }

AUTO_ASSERT( 0 == verify_data_segs_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ))
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

/* advance boundary seg pointer as an iterator.  returns NULL when iterator
   reaches end of list */
np02_boundary_seg *np02_boundary::advance_boundary_seg_itr(
    np02_boundary_seg *bndry_seg_itr) const{
AUTO_ASSERT( 0 == verify_data_segs_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ))
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

int np02_boundary::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const{
int err_cnt = 0;

if( NULL != m_shp_alloc ){
    const np02_boundary * const sa_bdry = 
        m_shp_alloc->alloc_get_boundary_by_idx(m_alloc_idx);
    if( this != sa_bdry ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary: this=%x  != (shp_alloc=%x)->boundary_by_idx(i=%i)=%x\n",
            this, m_shp_alloc, m_alloc_idx, sa_bdry );
        }
    }

if( NULL != m_owner ){
    const np02_shp_alloc *ownr_sa = m_owner->get_shp_alloc();
    if( ownr_sa != m_shp_alloc ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary: this=%x  shp_alloc=%x  !=  (owner=%x)->shp_alloc=%x\n",
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
            "boundary: this=%x  !=  (boundary_seg=%x)->owner=%x\n",
            this, s, s->get_owner());
        }
    if( m_shp_alloc != s->get_shp_alloc() ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary: this=%x shp_alloc=%x  !=  "
            "(boundary_seg=%x)->shp_alloc=%x\n",
            this, m_shp_alloc, s, s->get_shp_alloc());
        }
    s = advance_boundary_seg_citr(s);
    }

if( NULL != m_rgn_bdry_prev ){
    if( this != m_rgn_bdry_prev->m_rgn_bdry_next ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary: this=%x  !=  (rgn_bdry_prev=%x)->rgn_bdry_next=%x\n",
            this, m_rgn_bdry_prev, m_rgn_bdry_prev->m_rgn_bdry_next );
        }
    if( m_owner != m_rgn_bdry_prev->m_owner ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary:(this=%x)->owner=%x  !=  (rgn_bdry_prev=%x)->owner=%x\n",
            this, m_owner, m_rgn_bdry_prev, m_rgn_bdry_prev->m_owner );
        }
    if( m_shp_alloc != m_rgn_bdry_prev->m_shp_alloc ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary:(this=%x)->shp_alloc=%x  !=  "
            "(rgn_bdry_prev=%x)->shp_alloc=%x\n",
            this, m_shp_alloc, m_rgn_bdry_prev, m_rgn_bdry_prev->m_shp_alloc );
        }
    }

if( NULL != m_rgn_bdry_next ){
    if( this != m_rgn_bdry_next->m_rgn_bdry_prev ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary: this=%x  !=  (rgn_bdry_next=%x)->rgn_bdry_prerv=%x\n",
            this, m_rgn_bdry_next, m_rgn_bdry_next->m_rgn_bdry_prev );
        }
    if( m_owner != m_rgn_bdry_next->m_owner ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary:(this=%x)->owner=%x  !=  (rgn_bdry_next=%x)->owner=%x\n",
            this, m_owner, m_rgn_bdry_next, m_rgn_bdry_next->m_owner );
        }
    if( m_shp_alloc != m_rgn_bdry_next->m_shp_alloc ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "boundary:(this=%x)->shp_alloc=%x  !=  "
            "(rgn_bdry_next=%x)->shp_alloc=%x\n",
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
            "boundary: this=%x  segs_head NULL  segs_tail=%x\n",
            this, m_segs_tail );
        }
    } 
else if ( NULL == m_segs_tail ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "boundary: this=%x  segs_head=%x  segs_tail NULL\n",
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
                    "boundary: this=%x  segs_head=segs_tail=%x  "
                    "prev NULL   next=%x\n",
                     this, m_segs_head, m_segs_head->get_next() );
                }
            } 
        else if(m_segs_head->get_next() == NULL){
                ++err_cnt;
                np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                    "boundary: this=%x  segs_head=segs_tail=%x  "
                    "prev=%x   next NULL\n",
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
                "boundary: this=%x  segs_head=%x != segs_tail=%x  "
                "segs_head->prev=%x\n",
                 this, m_segs_head, m_segs_tail, m_segs_head->get_prev() );
            }
        if(m_segs_head->get_next() == NULL ){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "boundary: this=%x  segs_head=%x != segs_tail=%x  "
                "segs_head->next NULL\n",
                 this, m_segs_head, m_segs_tail );
            }
        if(m_segs_tail->get_prev() == NULL ){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "boundary: this=%x  segs_head=%x != segs_tail=%x  "
                "segs_tail->prev NULL\n",
                 this, m_segs_head, m_segs_tail );
            }
        if(m_segs_tail->get_next() != NULL ){
            ++err_cnt;
            np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "boundary: this=%x  segs_head=%x != segs_tail=%x  "
                "segs_tail->next=%x\n",
                 this, m_segs_head, m_segs_tail, m_segs_tail->get_next() );
            }
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

void np02_boundary::write_dxf_file(const std::string& layer,
    const uint8_t& color, np02_dxf_file *dxf_file) const{
const np02_boundary_seg *s = m_segs_head;
while ( NULL != s ){
    s->write_dxf_file(layer, color, dxf_file);
    s = advance_boundary_seg_citr(s);
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
    AUTO_ASSERT( (NULL == n) != (m_boundaries_tail == b) );
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
AUTO_ASSERT( 0 == verify_data_boundaries_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ))
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
        }
    }
AUTO_ASSERT( 0 == verify_data_boundaries_head_tail(AA_ERR_BUF(),
    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR() ))
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

int np02_region::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const{
int err_cnt = 0;

if( NULL != m_shp_alloc ){
    const np02_region * sa_rgn =
        m_shp_alloc->alloc_get_region_by_idx(m_alloc_idx);
    if( this != sa_rgn ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "region: this=%x  !=  "
            "(shp_alloc=%x)->get_region_by_idx(i=%i)=%x\n",
            this, m_shp_alloc, m_alloc_idx, sa_rgn );
        }
    }

if( NULL != m_owner ){
    const np02_region * ownr_rgn = m_owner->get_region();
    if( this != ownr_rgn ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "region: this=%x  !=  (owner=%x)->region=%x\n",
            this, m_owner, ownr_rgn );
        }
    }

const np02_boundary *b = m_boundaries_head;
while ( NULL != b ) {
    err_cnt += b->verify_data( err_msg, err_msg_capacity, err_msg_pos );
    if( this != b->get_owner() ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "region: this=%x  !=  (boundary=%x)->owner=%x\n",
            this, b, b->get_owner() );
        }
    if( m_shp_alloc != b->get_shp_alloc() ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "region: (this=%x) shp_alloc=%x !=  (boundary=%x)->owner=%x\n",
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
            "region: this=%x  boundaries_head NULL  boundaries_tail=%x\n",
            this, m_boundaries_tail );
        }
    } 
else if ( NULL == m_boundaries_tail ){
    ++err_cnt;
    np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "region: this=%x  boundaries_head=%x  boundaries_tail NULL\n",
         this, m_boundaries_head );
    }
else{
    /* open (non-loop) sequence, size > 1 */
    if(m_boundaries_head->get_rgn_bdry_prev() != NULL ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "region: this=%x  boundaries_head=%x != boundaries_tail=%x  "
            "boundaries_head->rgn_bdry_prev=%x\n",
             this, m_boundaries_head, m_boundaries_tail,
             m_boundaries_head->get_rgn_bdry_prev() );
        }
    if(m_boundaries_head->get_rgn_bdry_next() == NULL ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "region: this=%x  boundaries_head=%x != boundaries_tail=%x  "
            "boundaries_head->rgn_bdry_next NULL\n",
             this, m_boundaries_head, m_boundaries_tail );
        }
    if(m_boundaries_tail->get_rgn_bdry_prev() == NULL ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "region: this=%x  boundaries_head=%x != boundaries_tail=%x  "
            "boundaries_tail->rgn_bdry_prev NULL\n",
             this, m_boundaries_head, m_boundaries_tail );
        }
    if(m_boundaries_tail->get_rgn_bdry_next() != NULL ){
        ++err_cnt;
        np02_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "region: this=%x  boundaries_head=%x != boundaries_tail=%x  "
            "boundaries_tail->rgn_bdry_next=%x\n",
             this, m_boundaries_head, m_boundaries_tail,
             m_boundaries_tail->get_rgn_bdry_next() );
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
        np02_bmp_file *bmp_file) const{
/* TODO: implement */
/* 
get locator grid 
get bounding box 
get bounding box indices in bmp coordinates (include at least one extra pixel on all sides )
push bb indices in stack 
while stack not empty 
   pop box indices from stack
   if box too big, split in half, put each half in stack
   else
       find box center
       find nearest region shape to center of box
       if near point of region shape is outside of box, or box size = 1x1,
           color entire box if point is on IN side of shape
       else
           split box in half, put each half in stack
*/

}

void np02_region::write_dxf_file(const std::string& layer,
        const uint8_t& color, np02_dxf_file *dxf_file) const{
const np02_boundary *b = m_boundaries_head;
while ( NULL != b ) {
    b->write_dxf_file(layer, color, dxf_file);
    b = b->get_rgn_bdry_next();
    }
}

} /* namespace np02 */
