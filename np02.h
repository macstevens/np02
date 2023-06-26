/* np02.h  Newport 02 Library

Copyright (c) 2023 Mac Stevens <stevensm@earthlink.net> <www.macstevens.net>

Permission to use, copy, modify, and distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

Reference: https://opensource.org/licenses/ISC
*/

#ifndef NP02_H
#define NP02_H

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include "np02_bmp.h"
#include "np02_dxf.h"

namespace np02 {


template <class T, uint32_t initial_capacity>
class np02_small_vec {
public:
    typedef typename std::vector<T>::size_type size_type;
private:
    T m_vec[initial_capacity];
    size_type m_sz;
    std::vector<T> *m_std_vec;
public:
    np02_small_vec():m_sz(0),m_std_vec(NULL){}
   ~np02_small_vec(){clear();}
    size_type size() const{
        return (NULL == m_std_vec) ? m_sz : m_std_vec->size();}
    T& at(size_type i){
        if(NULL == m_std_vec)
            { return m_vec[i]; }
        else
            { return m_std_vec->at(i); }
        }
    const T& at(size_type i) const{
        if(NULL == m_std_vec)
            { return m_vec[i]; }
        else
            { return m_std_vec->at(i); }
        }
    void push_back(const T& val){reserve(size()+1);
        if(NULL == m_std_vec){m_vec[m_sz]=val;++m_sz;}
        else{m_std_vec->push_back(val);}}
    void reserve(size_type c){
        if(NULL == m_std_vec){
            if(c>initial_capacity){
                m_std_vec=new std::vector<T>();
                m_std_vec->reserve(c);
                for(size_type i=0;i<m_sz;++i) {m_std_vec->push_back(m_vec[i]);}
                }
            }
        else{m_std_vec->reserve(c);}
        }
    void erase_at(size_type i){
        if(i < size()){
            if(NULL == m_std_vec){
                for(size_type k=i;(k+1)<m_sz;++k){m_vec[k]=m_vec[k+1];}
                --m_sz;}
            else{m_std_vec->erase(m_std_vec->begin()+i);}
            }
        }
    void clear(){if(NULL!= m_std_vec){delete m_std_vec;} m_sz=0;}
};

/* trim whitespace off ends */
inline void np02_trim_str(std::string *s){
size_t first_non_space_pos, last_non_space_pos;
first_non_space_pos = s->find_first_not_of(" \t\v\n\r");
last_non_space_pos = s->find_last_not_of(" \t\v\n\r");
if(first_non_space_pos == std::string::npos){
    s->erase();
    }
else{
    if(first_non_space_pos > 0){ 
        size_t len = (last_non_space_pos + 1) - first_non_space_pos;
        s->assign(s->substr(first_non_space_pos, len));
        }
    else if((last_non_space_pos+1) != s->length()){
        s->resize(last_non_space_pos+1);
        }
    }
}


/* np02_shape.cpp */
class np02_xy;
class np02_shape;
class np02_circle;
class np02_arc;
class np02_line_seg;
class np02_rect;
class np02_polygon;
class np02_loc_grid_dim;
class np02_loc_grid_node;
class np02_loc_grid;
class np02_shp_alloc;
class np02_shape_test;



/* STL containers*/
typedef std::pair<uint16_t,uint16_t> np02_uint16_pair;
typedef std::vector<np02_uint16_pair> np02_uint16_pair_vec;
typedef np02_uint16_pair_vec::iterator np02_uint16_pair_vec_itr;
typedef np02_uint16_pair_vec::const_iterator np02_uint16_pair_vec_citr;


size_t np02_str_to_size_t( const std::string& str,
    int *error_count, std::string *err_msg );
double np02_str_to_double( const std::string& str,
    int *error_count, std::string *err_msg );


class np02_xy{
private:
    double m_x,m_y;
public:
    np02_xy():m_x(0.0),m_y(0.0){}
    np02_xy(const double& x, const double& y):m_x(x),m_y(y){}
    np02_xy(const np02_xy& master):m_x(master.m_x),m_y(master.m_y){}
   ~np02_xy(){}
    np02_xy& operator=(const np02_xy& master){
        m_x=master.m_x; m_y=master.m_y; return *this;}
    bool operator==(const np02_xy& other)const{
        return((other.m_x == m_x) && (other.m_y == m_y));}
    bool operator!=(const np02_xy& other)const{
        return((other.m_x != m_x) || (other.m_y != m_y));}
    bool operator<(const np02_xy& other)const{
        return((other.m_x<m_x) || ((other.m_x==m_x) && (other.m_y<m_y)));}
    const double& get_x() const { return m_x; }
    const double& get_y() const { return m_y; }
    void set_x(const double& x){m_x = x;}
    void set_y(const double& y){m_y = y;}
    void set_xy(const double& x, const double& y){m_x = x; m_y = y;}
    double get_len_sq() const{
        const double len_sq = ( m_x * m_x ) + ( m_y * m_y );
        return len_sq; }
    double get_len() const{
        const double len = sqrt( get_len_sq() );
        return len; }
    double get_dsq_to(const np02_xy& xy) const{
        const double dx = m_x - xy.get_x();
        const double dy = m_y - xy.get_y();
        return ((dx*dx) + (dy*dy)); }
    double get_distance_to(const np02_xy& xy) const{
        const double dx = m_x - xy.get_x();
        const double dy = m_y - xy.get_y();
        return sqrt((dx*dx) + (dy*dy)); }
    np02_xy get_unit_vector_to(const np02_xy& b) const{
        const double dx = b.get_x() - m_x;
        const double dy = b.get_y() - m_y;
        np02_xy u(1.0, 0.0);
        const double d_sq = (dx*dx) + (dy*dy); 
        if( d_sq > 0.0 ){
            const double d = sqrt(d_sq);
            u.set_x(dx/d);
            u.set_y(dy/d);
            }
        return u;
        }
    double cross(const np02_xy& xy) const{
        const double c = (m_x * xy.get_y()) - (m_y * xy.get_x());
        return c; }
    double dot(const np02_xy& xy) const{
        const double d = (m_x * xy.get_x()) + (m_y * xy.get_y());
        return d; }
    uint64_t hash( const uint64_t& h_in=0 ) const;
};

typedef std::vector<np02_xy> np02_xy_vec;
typedef np02_xy_vec::const_iterator np02_xy_vec_citr;
typedef np02_xy_vec::iterator np02_xy_vec_itr;

typedef uint32_t shape_idx_type;
typedef uint64_t shp_owner_idx_type;
typedef uint64_t lyr_idx_type;

#define NP02_SHAPE_MAX_IDX (std::numeric_limits<shape_idx_type>::max()/4)
#define NP02_SHP_OWNER_INVALID_IDX (std::numeric_limits<shp_owner_idx_type>::max())
#define NP02_LYR_INVALID_IDX (std::numeric_limits<lyr_idx_type>::max())



enum np02_answer_quality{
    NP02_ANSWER_QUALITY_SMALL_ERROR = 0,
    NP02_ANSWER_QUALITY_MEDIUM_ERROR = 1,
    NP02_ANSWER_QUALITY_LARGE_ERROR = 2,
    NP02_ANSWER_QUALITY_COUNT
};

/* helper struct for near point calculation*/
struct np02_dist_from_xy_xy{
public:
    np02_answer_quality m_answer_quality;
    double m_distance_from;
    bool m_near_xy_defined;
    bool m_other_near_xy_defined;
    np02_xy m_near_xy;
    np02_xy m_other_near_xy;
public:
    np02_dist_from_xy_xy(): m_answer_quality(NP02_ANSWER_QUALITY_COUNT),
        m_distance_from(0.0), m_near_xy_defined(false),
        m_other_near_xy_defined(false), m_near_xy(0.0,0.0),
        m_other_near_xy(0.0,0.0){}
    np02_dist_from_xy_xy(const np02_dist_from_xy_xy& master): 
        m_answer_quality(master.m_answer_quality),
        m_distance_from(master.m_distance_from), 
        m_near_xy_defined(master.m_near_xy_defined), 
        m_other_near_xy_defined(master.m_other_near_xy_defined),
        m_near_xy(master.m_near_xy), 
        m_other_near_xy(master.m_other_near_xy){}
   ~np02_dist_from_xy_xy(){}
    np02_dist_from_xy_xy& operator=(const np02_dist_from_xy_xy& master){
        m_answer_quality = master.m_answer_quality;
        m_distance_from = master.m_distance_from;
        m_near_xy_defined = master.m_near_xy_defined; 
        m_other_near_xy_defined = master.m_other_near_xy_defined;
        m_near_xy = master.m_near_xy; 
        m_other_near_xy = master.m_other_near_xy;
        return *this;}
        /* result < 0  ==> this is a better (more correct) result than other 
           result > 0  ==> this is a worse (less correct)result than other */
    int compare(const np02_dist_from_xy_xy& other) const{
        int result = static_cast<int>(m_answer_quality) -
            static_cast<int>(other.m_answer_quality);
        if(0 == result){
            if( m_distance_from < other.m_distance_from){ result = -1; }
            else if( other.m_distance_from < m_distance_from){ result = 1; }  }
        if(0 == result){ result = static_cast<int>(other.m_near_xy_defined) -
            static_cast<int>(m_near_xy_defined); }
        if(0 == result){ 
            result = static_cast<int>(other.m_other_near_xy_defined) -
            static_cast<int>(m_other_near_xy_defined); }
        if(0 == result){
            if( m_near_xy < other.m_near_xy){ result = -1; }
            else if( other.m_near_xy < m_near_xy){ result = 1; }  }
        if(0 == result){
            if( m_other_near_xy < other.m_other_near_xy){ result = -1; }
            else if( other.m_other_near_xy < m_other_near_xy){ result = 1; }  }
        return result; }
    bool operator==(const np02_dist_from_xy_xy& other) const{
        const int compare_result = compare(other);
        return (0 == compare_result) ? true : false; }
    bool operator<(const np02_dist_from_xy_xy& other) const{
        const int compare_result = compare(other);
        return (compare_result < 0) ? true : false; }
    void swap_xy(){
        const bool temp_defined = m_near_xy_defined;
        m_near_xy_defined = m_other_near_xy_defined;
        m_other_near_xy_defined = temp_defined;
        const np02_xy temp_xy = m_near_xy;
        m_near_xy = m_other_near_xy;
        m_other_near_xy = temp_xy; }
};

enum np02_rect_pt_idx{
    NP02_RECT_PT_IDX_P00 = 0,
    NP02_RECT_PT_IDX_P01 = 1,
    NP02_RECT_PT_IDX_P10 = 2,
    NP02_RECT_PT_IDX_P11 = 3,
    NP02_RECT_PT_IDX_CTR = 4,
    NP02_RECT_PT_IDX_COUNT
};

enum np02_shape_type{
    NP02_SHAPE_TYPE_SHAPE,
    NP02_SHAPE_TYPE_CIRCLE,
    NP02_SHAPE_TYPE_ARC,
    NP02_SHAPE_TYPE_LINE_SEG,
    NP02_SHAPE_TYPE_RECT,
    NP02_SHAPE_TYPE_POLYGON,
    NP02_SHAPE_TYPE_COUNT
};

class np02_shape{
public:
    static const double m_pi;
    static const double m_degrees_per_radian; 
    static const double m_little_ratio;  
    static const double m_little_ratio_sq; 
    static const double m_small_ratio;
    static const double m_small_ratio_sq;
    static const double m_tiny_ratio;
    static const double m_teensy_ratio;
    static const double m_teensy_ratio_sq;
    static const double m_googol;
public:
    typedef std::pair<np02_shape_type,uint32_t> type_idx_pair;
    typedef std::pair<type_idx_pair,np02_shape*> type_idx_shp;
    typedef std::vector<type_idx_shp> type_idx_shp_vec;
    typedef type_idx_shp_vec::iterator type_idx_shp_vec_itr;
    typedef type_idx_shp_vec::const_iterator type_idx_shp_vec_citr;
private:
    np02_shape_type m_shape_type; /* type+idx uniquely identifies shape */
    shape_idx_type m_shape_idx;
    shp_owner_idx_type m_shp_owner_idx;
    np02_shp_alloc *m_shp_alloc;
    np02_loc_grid_node *m_head_loc_grid_node;
protected:
    void set_shape_type(const np02_shape_type& shape_type){
        m_shape_type = shape_type; }
    double get_small_distance() const;
public:
    static void cos_sin_rot_deg( const double& rot_deg,
        double *cos_rot, double *sin_rot );
public:
    np02_shape();
    virtual void destruct();
    virtual ~np02_shape();
    np02_shape_type get_shape_type() const{return m_shape_type;}
    void set_shape_idx(const shape_idx_type& shape_idx){m_shape_idx=shape_idx;}
    shape_idx_type get_shape_idx() const{return m_shape_idx;}
    shp_owner_idx_type get_shp_owner_index() const{ return m_shp_owner_idx; }
    void set_shp_owner_idx(const shp_owner_idx_type& i){m_shp_owner_idx=i;}
    void set_shp_alloc(np02_shp_alloc *shp_alloc){ m_shp_alloc = shp_alloc; }
    np02_shp_alloc *get_shp_alloc() const { return m_shp_alloc; }
    lyr_idx_type get_lyr_idx() const;
    void set_head_loc_grid_node(np02_loc_grid_node *n){
        m_head_loc_grid_node = n; }
    np02_loc_grid_node *get_head_loc_grid_node() const{
        return m_head_loc_grid_node;}
    np02_loc_grid *get_loc_grid() const;
    void get_local_shapes(type_idx_shp_vec *local_shapes) const;

    virtual void get_bb(np02_xy *xy_min, np02_xy *xy_max) const=0;
    virtual void get_loc_grid_indices_for_init(
        const np02_loc_grid_dim& loc_grid_dim,
        const double& extra_search_d, np02_uint16_pair_vec *index_vec)
        const=0;
    virtual double get_distance_from_xy(const np02_xy& xy,
        np02_xy *near_xy=NULL) const=0;
    virtual double get_distance_from_line_seg_ab(const np02_xy& xy_a,
        const np02_xy& xy_b, np02_xy *near_xy=NULL,
        np02_xy *other_near_xy=NULL) const=0;
    virtual double get_distance_from_circle(const np02_circle *c,
        np02_xy *near_xy=NULL, np02_xy *circle_near_xy=NULL) const=0;
    virtual double get_distance_from_arc(const np02_arc *a,
        np02_xy *near_xy=NULL, np02_xy *arc_near_xy=NULL) const=0;
    virtual double get_distance_from_line_seg(const np02_line_seg *n,
        np02_xy *near_xy=NULL, np02_xy *line_seg_near_xy=NULL)const=0;
    virtual double get_distance_from_rect(const np02_rect *r,
        np02_xy *near_xy=NULL, np02_xy *rect_near_xy=NULL) const=0;
    virtual double get_distance_from_polygon(const np02_polygon *p,
        np02_xy *near_xy=NULL, np02_xy *polygon_near_xy=NULL)const=0;
    virtual double get_distance_from_shape(const np02_shape *s,
        np02_xy *near_xy=NULL, np02_xy *other_near_xy=NULL) const=0;
    double get_distance_from_shape_double_check(const np02_shape *s,
        double *err_estimate, np02_xy *near_xy=NULL,
        np02_xy *other_near_xy=NULL) const;
    virtual void translate(const np02_xy& dxy);
    virtual void translate_no_loc_grid(const np02_xy& dxy) = 0;
    virtual void rotate(const np02_xy& rot_ctr, const double& rot_deg);
    virtual void rotate_no_loc_grid(const np02_xy& rot_ctr,
        const double& rot_deg) = 0;

    virtual uint64_t hash( const uint64_t& h_in=0 ) const;
    virtual int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    static int verify_distance_from_shape_result(
        const np02_shape *shape_a,
        const np02_shape *shape_b, const np02_xy& xy_near_a,
        const np02_xy& xy_near_b, const double& distance_from,
        char *err_msg, const size_t err_msg_capacity, size_t *err_msg_pos );
    /* bounding box overlap 
           AAAAAAAAAA                             AAAAAAAAAA        
           min    max                             min    max        
                                                     
                   BBBBBBBBBB              BBBBBBBBBB
                   min    max              min    max
    */
    static bool is_bb_overlap(
        const np02_xy& bb_a_xy_min, const np02_xy& bb_a_xy_max,
        const np02_xy& bb_b_xy_min, const np02_xy& bb_b_xy_max ){
        const bool bb_overlap = ((bb_a_xy_min.get_x() <= bb_b_xy_max.get_x())&&
                                 (bb_b_xy_min.get_x() <= bb_a_xy_max.get_x())&&
                                 (bb_a_xy_min.get_y() <= bb_b_xy_max.get_y())&&
                                 (bb_b_xy_min.get_y() <= bb_a_xy_max.get_y()));
        return bb_overlap; }
    static double get_bb_sep_dist_manhattan(
        const np02_xy& bb_a_xy_min, const np02_xy& bb_a_xy_max,
        const np02_xy& bb_b_xy_min, const np02_xy& bb_b_xy_max );
    virtual std::ostream& ostream_output(std::ostream& os) const;
    virtual void write_bmp_file(const np02_xy& xy_min,
        const double& pixel_num, const np02_bmp_color& color,
        np02_bmp_file *bmp_file) const;
    virtual void write_dxf_file(const std::string& layer,
        const uint8_t& color, np02_dxf_file *dxf_file) const;
};

typedef std::vector<np02_shape *> np02_shape_vec;
typedef np02_shape_vec::iterator np02_shape_vec_itr;
typedef np02_shape_vec::const_iterator np02_shape_vec_citr;

uint64_t np02_shape_vec_hash( const np02_shape_vec& v, 
    const uint64_t& h_in = 0 );


class np02_circle: public np02_shape{
private:
    np02_xy m_ctr;
    double m_radius;
public:
    np02_circle();
    virtual ~np02_circle();
    void set_free_chain_next(np02_circle *n){
        set_head_loc_grid_node(reinterpret_cast<np02_loc_grid_node *>(n));}
    np02_circle *get_free_chain_next() const{ return
        reinterpret_cast<np02_circle *>(get_head_loc_grid_node()); }
    void init(const np02_xy& ctr, const double& radius){
        m_ctr=ctr; m_radius=radius; }
    void set_ctr(const np02_xy& ctr){m_ctr=ctr;}
    const np02_xy& get_ctr() const{return m_ctr;}
    void set_radius(const double& radius){assert(radius>=0.0);m_radius=radius;}
    const double& get_radius() const{return m_radius;}
    virtual void get_bb(np02_xy *xy_min, np02_xy *xy_max) const;
    virtual void get_loc_grid_indices_for_init(
        const np02_loc_grid_dim& loc_grid_dim,
        const double& extra_search_d, np02_uint16_pair_vec *index_vec)
        const;
    virtual double get_distance_from_xy(const np02_xy& xy,
        np02_xy *near_xy=NULL) const;
    virtual double get_distance_from_line_seg_ab(const np02_xy& xy_a,
        const np02_xy& xy_b, np02_xy *near_xy=NULL,
        np02_xy *other_near_xy=NULL) const;
    virtual double get_distance_from_circle(const np02_circle *c,
        np02_xy *near_xy=NULL, np02_xy *circle_near_xy=NULL) const;
    virtual double get_distance_from_arc(const np02_arc *a,
        np02_xy *near_xy=NULL, np02_xy *arc_near_xy=NULL) const;
    virtual double get_distance_from_line_seg(const np02_line_seg *n,
        np02_xy *near_xy=NULL, np02_xy *line_seg_near_xy=NULL)const;
    virtual double get_distance_from_rect(const np02_rect *r,
        np02_xy *near_xy=NULL, np02_xy *rect_near_xy=NULL) const;
    virtual double get_distance_from_polygon(const np02_polygon *p,
        np02_xy *near_xy=NULL, np02_xy *polygon_near_xy=NULL)const;
    virtual double get_distance_from_shape(const np02_shape *s,
        np02_xy *near_xy=NULL, np02_xy *other_near_xy=NULL) const;
    virtual void translate_no_loc_grid(const np02_xy& dxy);
    virtual void rotate_no_loc_grid(const np02_xy& rot_ctr,
        const double& rot_deg);

    virtual uint64_t hash( const uint64_t& h_in=0 ) const;
    virtual int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    virtual std::ostream& ostream_output(std::ostream& os) const;
    virtual void write_bmp_file(const np02_xy& xy_min,
        const double& pixel_num, const np02_bmp_color& color,
        np02_bmp_file *bmp_file) const;
    virtual void write_dxf_file(const std::string& layer,
        const uint8_t& color, np02_dxf_file *dxf_file) const;

    size_t circle_circle_intersect(const np02_circle *c,
        np02_xy *xy_intsct_0=NULL, np02_xy *xy_intsct_1=NULL) const;
    static size_t circle_circle_intsct(const np02_xy& ctr_a,
        const double& r_a, const np02_xy& ctr_b, const double& r_b,
        np02_xy *xy_intsct_0=NULL, np02_xy *xy_intsct_1=NULL);
    static int verify_data_circle_circle_intsct_result( const np02_xy& ctr_a,
        const double& r_a, const np02_xy& ctr_b, const double& r_b,
        np02_xy *xy_intsct_0, np02_xy *xy_intsct_1, const size_t& intsct_count,
        char *err_msg, const size_t err_msg_capacity, size_t *err_msg_pos );
};


/*
                 \fwd0
                   \
                     \
                ***    \
           *           * \
        *                 +P0
      *                  .        direction = CCW
     +P1               .
    /     .           .start_angle
   /  end_angle.    .
  /               +ctr
fwd1



*/
class np02_arc: public np02_shape{
public:
    struct init_params{
       np02_xy m_ctr;
       double m_radius; 
       double m_start_angle_deg; 
       double m_end_angle_deg;
       double m_width;
    };
    struct init3pt_params{
       np02_xy m_pa;
       np02_xy m_pb;
       np02_xy m_pc;
       double m_width;
       double m_max_radius;
    };
    struct init3pt_aux_params{
       np02_xy m_p_0;
       np02_xy m_p_1;
       bool m_is_straight;
    };
private:
    /* init data */
    np02_xy m_ctr;
    double m_radius;
    double m_start_angle_deg;
    double m_end_angle_deg;
    double m_width;

    /* derived values */
    np02_xy m_p_0, m_p_1;
    np02_xy m_fwd_0, m_fwd_1;
    double m_fwd_dot_0, m_fwd_dot_1;
    np02_xy m_bb_xy_min, m_bb_xy_max;
public:
    np02_arc();
    virtual ~np02_arc();
    void set_free_chain_next(np02_arc *n){
        set_head_loc_grid_node(reinterpret_cast<np02_loc_grid_node *>(n));}
    np02_arc *get_free_chain_next() const{ return
        reinterpret_cast<np02_arc *>(get_head_loc_grid_node()); }
    void init(const init_params& prm);
    void init_force_p01(const init_params& prm, 
        const init3pt_aux_params& aux_prm );
    static void init3pt_to_init_params(const init3pt_params& init3pt_prm,
        init_params *init_prm, init3pt_aux_params *aux_params ); 
    const np02_xy& get_ctr() const{return m_ctr;}
    const double& get_radius() const{return m_radius;}
    const double& get_start_angle_deg() const{return m_start_angle_deg;}
    const double& get_end_angle_deg() const{return m_end_angle_deg;}
    const double& get_width() const{return m_width;}
    const np02_xy& get_p_0() const{return m_p_0;}
    const np02_xy& get_p_1() const{return m_p_1;}
    const np02_xy& get_fwd_0() const{return m_fwd_0;}
    const np02_xy& get_fwd_1() const{return m_fwd_1;}
    const double& get_fwd_dot_0() const{return m_fwd_dot_0;}
    const double& get_fwd_dot_1() const{return m_fwd_dot_1;}
    bool is_less_than_half_circle() const{
        const double angle_range = m_end_angle_deg-m_start_angle_deg; 
        return ( angle_range < 180.0 ) ? true : false; }
    /* result > 0 => p is inside zone 0 */
    double distance_in_zone_0(const np02_xy& p) const{
        const double fwd_0_dot_p = m_fwd_0.dot(p);
        const double d_in_z_0 = m_fwd_dot_0 - fwd_0_dot_p;
        return d_in_z_0; }
    /* result > 0 => p is inside zone 1 */
    double distance_in_zone_1(const np02_xy& p) const{
        const double fwd_1_dot_p = m_fwd_1.dot(p);
        const double d_in_z_1 = fwd_1_dot_p - m_fwd_dot_1;
        return d_in_z_1; }
    bool is_in_zone_0(const np02_xy& p) const{
        const double fwd_0_dot_p = m_fwd_0.dot(p);
        return ( fwd_0_dot_p < m_fwd_dot_0 ) ? true : false; }
    bool is_in_zone_1(const np02_xy& p) const{
        const double fwd_1_dot_p = m_fwd_1.dot(p);
        return ( fwd_1_dot_p > m_fwd_dot_1 ) ? true : false; }
    bool is_in_perp_zone(const np02_xy& p) const;
    double distance_in_perp_zone(const np02_xy& p) const;
    virtual void get_bb(np02_xy *xy_min, np02_xy *xy_max) const;
    virtual void get_loc_grid_indices_for_init(
        const np02_loc_grid_dim& loc_grid_dim,
        const double& extra_search_d, np02_uint16_pair_vec *index_vec)
        const;
    virtual double get_distance_from_xy(const np02_xy& xy,
        np02_xy *near_xy=NULL) const;
    virtual double get_distance_from_line_seg_ab(const np02_xy& xy_a,
        const np02_xy& xy_b, np02_xy *near_xy=NULL,
        np02_xy *other_near_xy=NULL) const;
    virtual double get_distance_from_circle(const np02_circle *c,
        np02_xy *near_xy=NULL, np02_xy *circle_near_xy=NULL) const;
    virtual double get_distance_from_arc(const np02_arc *a,
        np02_xy *near_xy=NULL, np02_xy *arc_near_xy=NULL) const;
    virtual double get_distance_from_line_seg(const np02_line_seg *n,
        np02_xy *near_xy=NULL, np02_xy *line_seg_near_xy=NULL)const;
    virtual double get_distance_from_rect(const np02_rect *r,
        np02_xy *near_xy=NULL, np02_xy *rect_near_xy=NULL) const;
    virtual double get_distance_from_polygon(const np02_polygon *p,
        np02_xy *near_xy=NULL, np02_xy *polygon_near_xy=NULL)const;
    virtual double get_distance_from_shape(const np02_shape *s,
        np02_xy *near_xy=NULL, np02_xy *other_near_xy=NULL) const;
    virtual void translate_no_loc_grid(const np02_xy& dxy);
    virtual void rotate_no_loc_grid(const np02_xy& rot_ctr,
        const double& rot_deg);

    virtual uint64_t hash( const uint64_t& h_in=0 ) const;
    virtual int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    virtual std::ostream& ostream_output(std::ostream& os) const;
    virtual void write_bmp_file(const np02_xy& xy_min,
        const double& pixel_num, const np02_bmp_color& color,
        np02_bmp_file *bmp_file) const;
    virtual void write_dxf_file(const std::string& layer,
        const uint8_t& color, np02_dxf_file *dxf_file) const;

    size_t arc_circle_centerline_intersect(const np02_circle *c,
        np02_xy *xy_intsct_0=NULL, np02_xy *xy_intsct_1=NULL) const;
    int verify_data_arc_circle_centerline_intersect_result(
        const np02_circle *c, np02_xy *xy_intsct_0, np02_xy *xy_intsct_1,
        const size_t& intsct_count, char *err_msg,
        const size_t err_msg_capacity, size_t *err_msg_pos ) const;
    size_t arc_arc_centerline_intersect(const np02_arc *b,
        np02_xy *xy_intsct_0=NULL, np02_xy *xy_intsct_1=NULL) const;
    int verify_data_arc_arc_centerline_intersect_result( const np02_arc *a,
        np02_xy *xy_intsct_0, np02_xy *xy_intsct_1, const size_t& intsct_count,
        char *err_msg, const size_t err_msg_capacity,
        size_t *err_msg_pos ) const;
    size_t arc_seg_centerline_intersect(const np02_line_seg *n,
        np02_xy *xy_intsct_0=NULL, np02_xy *xy_intsct_1=NULL) const;
    int verify_data_arc_seg_centerline_intersect_result(
        const np02_line_seg *n, np02_xy *xy_intsct_0, np02_xy *xy_intsct_1,
        const size_t& intsct_count, char *err_msg,
        const size_t err_msg_capacity, size_t *err_msg_pos ) const;
private:
    void init_bb();
    double get_distance_from_xy_hw(const np02_xy& xy,
        const double& hw, np02_xy *near_xy=NULL) const;
    void get_distance_from_arc_zone_0(const np02_arc *a,
        np02_dist_from_xy_xy *result) const;
    void get_distance_from_arc_zone_1(const np02_arc *a,
        np02_dist_from_xy_xy *result) const;
    void get_distance_from_arc_far_perp_zone(const np02_arc *a,
        np02_dist_from_xy_xy *result) const;
    void get_distance_from_arc_centerline_intersect(const np02_arc *a,
        np02_dist_from_xy_xy *result_0, np02_dist_from_xy_xy *result_1) const;
    void get_distance_from_shape_zone_01(const np02_shape *shp,
        const bool& is_from_p0, np02_dist_from_xy_xy *result) const;
    void get_distance_from_shape_far_perp_zone(const np02_shape *shp,
        np02_dist_from_xy_xy *result) const;
    void get_distance_from_line_seg_far_perp_zone(const np02_line_seg *n,
        np02_dist_from_xy_xy *result) const;
    void get_distance_from_rect_pt(const np02_rect *rect,
        const np02_rect_pt_idx& rect_pt_idx, 
        np02_dist_from_xy_xy *result) const;
    void get_distance_from_line_seg_p01(const np02_line_seg *n,
        const bool& is_from_seg_p0, np02_dist_from_xy_xy *result) const;
    void get_distance_from_line_seg_centerline_intersect(const np02_line_seg *n,
        np02_dist_from_xy_xy *result_0, np02_dist_from_xy_xy *result_1) const;
};


/* rectangle 
                             fwd_dot_ctr
                             ---->
  fwd_cross_1 ^           p01             p11
              |           *---------------*
              |           |               |                               fwd 
                          |      ctr      |                             /         
      fwd_cross_ctr^     h|       +-------|--------->fwd              /           
                   |      |               |                         /             
                   |      |               |                       /   rot_deg     
  fwd_cross_0 ^           *---------------*                     +- - - - - - -  
              |           p00     w       p10
              |      
                     ---->            ---->
                 fwd_dot_0        fwd_dot_1
*/
class np02_rect: public np02_shape{
private:
    /* init data */
    np02_xy m_ctr;
    double m_w, m_h; /* width, height */
    double m_rot_deg; /* rotation, in degrees */

    /* derived values */
    np02_xy m_fwd;
    double m_fwd_dot_ctr, m_fwd_cross_ctr;
    np02_xy m_p00,m_p01,m_p10,m_p11; /* corners */
public:
    np02_rect();
    virtual ~np02_rect();
    void set_free_chain_next(np02_rect *n){
        set_head_loc_grid_node(reinterpret_cast<np02_loc_grid_node *>(n));}
    np02_rect *get_free_chain_next() const{ return
        reinterpret_cast<np02_rect *>(get_head_loc_grid_node()); }
    void init(const np02_xy& ctr, const double& w, const double& h,
        const double& rot_deg);

    const np02_xy& get_ctr() const {return m_ctr;}
    const double& get_w() const { return m_w; }
    const double& get_h() const { return m_h; }
    const double& get_rot_deg() const { return m_rot_deg; }

    /* derived values */
    const np02_xy& get_fwd() const { return m_fwd; }
    const double& get_fwd_dot_ctr() const { return m_fwd_dot_ctr; }
    const double get_fwd_dot_0() const { return (m_fwd_dot_ctr-(m_w/2.0)); }
    const double get_fwd_dot_1() const { return (m_fwd_dot_ctr+(m_w/2.0)); }
    const double& get_fwd_cross_ctr() const { return m_fwd_cross_ctr; }
    const double get_fwd_cross_0() const { return (m_fwd_cross_ctr-(m_h/2.0));}
    const double get_fwd_cross_1() const { return (m_fwd_cross_ctr+(m_h/2.0));} 
    const np02_xy& get_p00() const { return m_p00; }
    const np02_xy& get_p01() const { return m_p01; }
    const np02_xy& get_p10() const { return m_p10; }
    const np02_xy& get_p11() const { return m_p11; }
    const np02_xy& get_pt_by_idx( const np02_rect_pt_idx& i ) const;

    virtual void get_bb(np02_xy *xy_min, np02_xy *xy_max) const;
    virtual void get_loc_grid_indices_for_init(
        const np02_loc_grid_dim& loc_grid_dim,
        const double& extra_search_d, np02_uint16_pair_vec *index_vec)
        const;
    virtual double get_distance_from_xy(const np02_xy& xy,
        np02_xy *near_xy=NULL) const;
    virtual double get_distance_from_line_seg_ab(const np02_xy& xy_a,
        const np02_xy& xy_b, np02_xy *near_xy=NULL,
        np02_xy *other_near_xy=NULL) const;
    virtual double get_distance_from_circle(const np02_circle *c,
        np02_xy *near_xy=NULL, np02_xy *circle_near_xy=NULL) const;
    virtual double get_distance_from_arc(const np02_arc *a,
        np02_xy *near_xy=NULL, np02_xy *arc_near_xy=NULL) const;
    virtual double get_distance_from_line_seg(const np02_line_seg *n,
        np02_xy *near_xy=NULL, np02_xy *line_seg_near_xy=NULL)const;
    virtual double get_distance_from_rect(const np02_rect *r,
        np02_xy *near_xy=NULL, np02_xy *rect_near_xy=NULL) const;
    virtual double get_distance_from_polygon(const np02_polygon *p,
        np02_xy *near_xy=NULL, np02_xy *polygon_near_xy=NULL)const;
    virtual double get_distance_from_shape(const np02_shape *s,
        np02_xy *near_xy=NULL, np02_xy *other_near_xy=NULL) const;
    virtual void translate_no_loc_grid(const np02_xy& dxy);
    virtual void rotate_no_loc_grid(const np02_xy& rot_ctr,
        const double& rot_deg);

    virtual uint64_t hash( const uint64_t& h_in=0 ) const;
    virtual int verify_data( char *err_msg, const size_t err_msg_capacity,
        size_t *err_msg_pos ) const;
    int verify_data_num( char *err_msg, const size_t err_msg_capacity,
        size_t *err_msg_pos ) const;
    int verify_data_rect_loc_grid( char *err_msg,
        const size_t err_msg_capacity, size_t *err_msg_pos ) const;
    virtual std::ostream& ostream_output(std::ostream& os) const;
    virtual void write_bmp_file(const np02_xy& xy_min,
        const double& pixel_num, const np02_bmp_color& color,
        np02_bmp_file *bmp_file) const;
    virtual void write_dxf_file(const std::string& layer,
        const uint8_t& color, np02_dxf_file *dxf_file) const;
};


class np02_line_seg: public np02_shape{
private:
    /* init data */
    np02_xy m_p_0, m_p_1;
    double m_width;

    /* derived values */
    np02_xy m_fwd;
    double m_fwd_dot_0, m_fwd_dot_1, m_fwd_cross_01;
public:
    np02_line_seg();
    virtual ~np02_line_seg();
    void set_free_chain_next(np02_line_seg *n){
        set_head_loc_grid_node(reinterpret_cast<np02_loc_grid_node *>(n));}
    np02_line_seg *get_free_chain_next() const{ return
        reinterpret_cast<np02_line_seg *>(get_head_loc_grid_node()); }
    void init(const np02_xy& p_0,const np02_xy& p_1,const double& width);
    const np02_xy& get_p_0() const{ return m_p_0; }
    const np02_xy& get_p_1() const{ return m_p_1; }
    const double& get_width() const{ return m_width; }
    double get_len() const{ return (m_fwd_dot_1 - m_fwd_dot_0); }
    const np02_xy& get_fwd() const{ return m_fwd; }
    const double& get_fwd_dot_0() const{ return m_fwd_dot_0; }
    const double& get_fwd_dot_1() const{ return m_fwd_dot_1; }
    const double& get_fwd_cross_01() const{ return m_fwd_cross_01; }
    virtual void get_bb(np02_xy *xy_min, np02_xy *xy_max) const;
    virtual void get_loc_grid_indices_for_init(
        const np02_loc_grid_dim& loc_grid_dim,
        const double& extra_search_d, np02_uint16_pair_vec *index_vec)
        const;
    virtual double get_distance_from_xy(const np02_xy& xy,
        np02_xy *near_xy=NULL) const;
    double get_distance_from_xy_hw(const np02_xy& xy,
        const double& hw, np02_xy *near_xy=NULL) const;
    virtual double get_distance_from_line_seg_ab(const np02_xy& xy_a,
        const np02_xy& xy_b, np02_xy *near_xy=NULL,
        np02_xy *other_near_xy=NULL) const;
    virtual double get_distance_from_circle(const np02_circle *c,
        np02_xy *near_xy=NULL, np02_xy *circle_near_xy=NULL) const;
    virtual double get_distance_from_arc(const np02_arc *a,
        np02_xy *near_xy=NULL, np02_xy *arc_near_xy=NULL) const;
    virtual double get_distance_from_line_seg(const np02_line_seg *n,
        np02_xy *near_xy=NULL, np02_xy *line_seg_near_xy=NULL)const;
    virtual double get_distance_from_line_seg_double_check(
        const np02_line_seg *n, double *err_estimate,
        np02_xy *near_xy=NULL, np02_xy *line_seg_near_xy=NULL)const;
    virtual double get_distance_from_rect(const np02_rect *r,
        np02_xy *near_xy=NULL, np02_xy *rect_near_xy=NULL) const;
    virtual double get_distance_from_polygon(const np02_polygon *p,
        np02_xy *near_xy=NULL, np02_xy *polygon_near_xy=NULL)const;
    virtual double get_distance_from_shape(const np02_shape *s,
        np02_xy *near_xy=NULL, np02_xy *other_near_xy=NULL) const;
    virtual void translate_no_loc_grid(const np02_xy& dxy);
    virtual void rotate_no_loc_grid(const np02_xy& rot_ctr,
        const double& rot_deg);

    virtual uint64_t hash( const uint64_t& h_in=0 ) const;
    virtual int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    int verify_data_num( char *err_msg, const size_t err_msg_capacity,
        size_t *err_msg_pos ) const;
    virtual std::ostream& ostream_output(std::ostream& os) const;
    virtual void write_bmp_file(const np02_xy& xy_min,
        const double& pixel_num, const np02_bmp_color& color,
        np02_bmp_file *bmp_file) const;
    virtual void write_dxf_file(const std::string& layer,
        const uint8_t& color, np02_dxf_file *dxf_file) const;
};


class np02_polygon: public np02_shape{
public:
    np02_xy_vec m_vertices;
public:
    np02_polygon();
    virtual ~np02_polygon();
    void set_free_chain_next(np02_polygon *n){
        set_head_loc_grid_node(reinterpret_cast<np02_loc_grid_node *>(n));}
    np02_polygon *get_free_chain_next() const{ return
        reinterpret_cast<np02_polygon *>(get_head_loc_grid_node()); }
    void init(const np02_xy_vec& vertices);
    virtual void get_bb(np02_xy *xy_min, np02_xy *xy_max) const;
    virtual void get_loc_grid_indices_for_init(
        const np02_loc_grid_dim& loc_grid_dim,
        const double& extra_search_d, np02_uint16_pair_vec *index_vec)
        const;
    virtual double get_distance_from_xy(const np02_xy& xy,
        np02_xy *near_xy=NULL) const;
    virtual double get_distance_from_line_seg_ab(const np02_xy& xy_a,
        const np02_xy& xy_b, np02_xy *near_xy=NULL,
        np02_xy *other_near_xy=NULL) const;
    virtual double get_distance_from_circle(const np02_circle *c,
        np02_xy *near_xy=NULL, np02_xy *circle_near_xy=NULL) const;
    virtual double get_distance_from_arc(const np02_arc *a,
        np02_xy *near_xy=NULL, np02_xy *arc_near_xy=NULL) const;
    virtual double get_distance_from_line_seg(const np02_line_seg *n,
        np02_xy *near_xy=NULL, np02_xy *line_seg_near_xy=NULL)const;
    virtual double get_distance_from_rect(const np02_rect *r,
        np02_xy *near_xy=NULL, np02_xy *rect_near_xy=NULL) const;
    virtual double get_distance_from_polygon(const np02_polygon *p,
        np02_xy *near_xy=NULL, np02_xy *polygon_near_xy=NULL)const;
    virtual double get_distance_from_shape(const np02_shape *s,
        np02_xy *near_xy=NULL, np02_xy *other_near_xy=NULL) const;
    virtual void translate_no_loc_grid(const np02_xy& dxy);
    virtual void rotate_no_loc_grid(const np02_xy& rot_ctr,
        const double& rot_deg);

    virtual uint64_t hash( const uint64_t& h_in=0 ) const;
    virtual int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    virtual std::ostream& ostream_output(std::ostream& os) const;
    virtual void write_bmp_file(const np02_xy& xy_min,
        const double& pixel_num, const np02_bmp_color& color,
        np02_bmp_file *bmp_file) const;
    virtual void write_dxf_file(const std::string& layer,
        const uint8_t& color, np02_dxf_file *dxf_file) const;
};


/* initialize from points */
struct np02_loc_grid_dim_init_params{
    size_t point_count; /* number of points */
    np02_xy bb_min_xy, bb_max_xy;  /* bounding box of points */
    double loc_grid_density; /* target density = points / grid square */
    size_t max_loc_grid_sq_count; /* target max number of grid squares */
};


/* locator grid dimensions 

           < - - - - - - - - - w - - - - - - - - - >
                      (number of squares)
           +-------+-------+-------+-------+-------+-   ^
           |       |       |       |       |       |
           |       |       |       |       |       |    |
           |       |       |       |       |       |
           +-------+-------+-------+-------+-------+-   |
           |       |       |       |       |       |     
           |       |       |       |       |       |    h (number of squares)
           |       |       |       |       |       |
       ^   +-------+-------+-------+-------+-------+-   |
       |   |       |       |       |       |       |
     sq_sz |       |       |       |       |       |    |
       |   |       |       |       |       |       |
       v   +-------+-------+-------+-------+-------+-   v
           (x_min, y_min)                  <-sq_sz->
*/
class np02_loc_grid_dim {
private:
    uint16_t m_w, m_h; /* width, height */
    double m_x_min, m_y_min; /* lower left corner */
    double m_sq_size; /* square size */
public:
    np02_loc_grid_dim(): m_w(0), m_h(0), m_x_min(0.0), m_y_min(0.0),
        m_sq_size(0.0){}
    np02_loc_grid_dim(const np02_loc_grid_dim& d): m_w(d.m_w), 
        m_h(d.m_h), m_x_min(d.m_x_min), m_y_min(d.m_y_min),
        m_sq_size(d.m_sq_size){}
   ~np02_loc_grid_dim() {}

    np02_loc_grid_dim& operator=(const np02_loc_grid_dim& d)
        { m_w=d.m_w; m_h=d.m_h; m_x_min=d.m_x_min; m_y_min=d.m_y_min;
        m_sq_size=d.m_sq_size; return *this; }
    bool operator==(const np02_loc_grid_dim& other) const {
        return ((m_w == other.m_w) && (m_h == other.m_h) &&
            (m_x_min == other.m_x_min) && (m_y_min == other.m_y_min) &&
            (m_sq_size == other.m_sq_size)); }
    bool operator!=(const np02_loc_grid_dim& other) const {
        return ((m_w != other.m_w) && (m_h != other.m_h) &&
            (m_x_min != other.m_x_min) && (m_y_min != other.m_y_min) &&
            (m_sq_size != other.m_sq_size)); }

    void reset(){ m_w=0; m_h=0; m_x_min=0.0; m_y_min=0.0; m_sq_size=0.0; }
    void init(const np02_loc_grid_dim_init_params *init_params);
    void set_w(const uint16_t& w){ m_w = w; }
    void set_h(const uint16_t& h){ m_h = h; }
    void set_x_min(const double& x_min){ m_x_min = x_min; }
    void set_y_min(const double& y_min){ m_y_min = y_min; }
    void set_sq_size(const double& sq_size){ m_sq_size = sq_size; }

    uint16_t get_w() const { return m_w; }
    uint16_t get_h() const { return m_h; }
    const double& get_x_min() const { return m_x_min; }
    const double& get_y_min() const { return m_y_min; }
    const double& get_sq_size() const { return m_sq_size; }

    uint16_t get_i(const double& x) const {
        uint16_t i = 0;
        if((m_sq_size>0.0) && (m_w > 0)){
            double dx = x - m_x_min;  double ii = dx / m_sq_size;
            i= (ii<1.0)? 0 : (ii>=(double)m_w)? m_w-1: (uint16_t)ii; }
        return i; }
    uint16_t get_j(const double& y) const {
        uint16_t j = 0;
        if((m_sq_size>0.0) && (m_h > 0)){
            double dy = y - m_y_min;  double jj = dy / m_sq_size;
            j= (jj<1.0)? 0 : (jj>=(double)m_h)? m_h-1: (uint16_t)jj; }
        return j; }
    void get_bb_indices( const np02_xy& xy_min, const np02_xy& xy_max,
        np02_uint16_pair *ij_min, np02_uint16_pair *ij_max) const;
    double get_sq_ctr_x( const uint16_t& i ) const{
        double ctr_x;
        if(m_w>0){
            double ii=static_cast<double>((i>=m_w)?(m_w-1):i);
            ctr_x = m_x_min + ( m_sq_size * ( ii + 0.5 ) ); }
        else{ ctr_x = m_x_min; }
        return ctr_x; }
    double get_sq_ctr_y( const uint16_t& j ) const{
        double ctr_y;
        if(m_h>0){
            double jj=static_cast<double>((j>=m_h)?(m_h-1):j);
            ctr_y = m_y_min + ( m_sq_size * ( jj + 0.5 ) ); }
        else{ ctr_y = m_y_min; }
        return ctr_y; }

    uint64_t hash( const uint64_t& h_in=0 ) const;
    int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    std::ostream& ostream_output(std::ostream& os) const;
    void write_bmp_file(const np02_xy& xy_min,
        const double& pixel_num, const np02_bmp_color& color,
        np02_bmp_file *bmp_file) const;
};

/*
  +---------------------------------------------------------------------------------------+
  |np02_shape                                                                             |
  |m_loc_grid_head_node+                                                                  |
  +--------------------|------------------------------------------------------------------+
                       |  ^                   ^                   ^                   ^
                       v  |                   |                   |                   |
                  +-------|-----+     +-------|-----+     +-------|-----+     +-------|-----+
                  |m_owner+     |     |m_owner+     |     |m_owner+     |     |m_owner+     |
                  |m_s_prev=NULL|<-----m_s_prev     |<-----m_s_prev     |<-----m_s_prev     |
                  |m_s_next---------->|m_s_next---------->|m_s_next---------->|m_s_next=NULL|
                  +-------------+     +-------------+     +-------------+     +-------------+




  +----------------------------------------------------------------------------------------+
  |np02_loc_grid                                                                           |
  |            +-----------------+-----------------+-----------------+-----------------+   |                                                   |
  |            |   i=0     j=0   |   i=0     j=1   |   i=1     j=0   |   i=1     j=1   |   |
  |            +-----------------+-----------------+-----------------+-----------------+   |                                                   |
  | m_loc_grid | 0    *          | 1    *          | 2    *          | 3    *          |   |
  |            +------|----------+------|----------+------|----------+------|----------+   |                                                                         |
  +-------------------|-----------------|-----------------|-----------------|--------------+
                      |                 |                 |                 |
                      v                 v                 v                 v
           +------------+    +------------+    +------------+    +------------+
           |i=0 j=0     |    |i=0 j=1     |    |i=1 j=0     |    |i=1 j=1     |  
           |m_owner=1   |    |m_owner=1   |    |m_owner=2   |    |m_owner=1   |  
           |m_prev=NULL |    |m_prev=NULL |    |m_prev=NULL |    |m_prev=NULL |
           |m_next-+    |    |m_next-+    |    |m_next-+    |    |m_next-+    |
           +-------|----+    +-------|----+    +-------|----+    +-------|----+
                   |  ^              |  ^              |  ^              |  ^
                   v  |              v  |              v  |              v  |
           +----------|-+    +----------|-+    +----------|-+    +----------|-+
           |i=0 j=0   | |    |i=0 j=1   | |    |i=1 j=0   | |    |i=1 j=1   | |  
           |m_owner=2 | |    |m_owner=2 | |    |m_owner=3 | |    |m_owner=3 | |  
           |m_prev----+ |    |m_prev----+ |    |m_prev----+ |    |m_prev----+ |
           |m_next=NULL |    |m_next-+    |    |m_next-+    |    |m_next-+    |
           +------------+    +-------|----+    +-------|----+    +-------|----+
                                     |  ^              |  ^              |  ^  
                                     v  |              v  |              v  |  
                             +----------|-+    +----------|-+    +----------|-+
                             |i=0 j=1   | |    |i=1 j=0   | |    |i=1 j=1   | |
                             |m_owner=3 | |    |m_owner=4 | |    |m_owner=5 | |
                             |m_prev----+ |    |m_prev----+ |    |m_prev----+ |
                             |m_next-+    |    |m_next=NULL |    |m_next=NULL |
                             +-------|----+    +------------+    +------------+
                                     |  ^     
                                     v  |     
                             +----------|-+   
                             |i=0 j=1   | |   
                             |m_owner=4 | |   
                             |m_prev----+ |   
                             |m_next=NULL |   
                             +------------+   
                 
*/
class np02_loc_grid_node{
private:
    size_t m_alloc_idx;
    np02_shape *m_owner;
    np02_loc_grid *m_loc_grid;
    uint16_t m_i, m_j;
    np02_loc_grid_node *m_prev, *m_next; /* same grid square */
    np02_loc_grid_node *m_s_prev, *m_s_next; /* same shape */
public:
    np02_loc_grid_node(): m_alloc_idx(0), m_owner(NULL), m_loc_grid(NULL),
        m_i(0), m_j(0), m_prev(NULL), m_next(NULL), m_s_prev(NULL),
        m_s_next(NULL) {}
    ~np02_loc_grid_node(){}

    void set_free_chain_next(np02_loc_grid_node *n){ m_next = n;}
    np02_loc_grid_node *get_free_chain_next() const{ return m_next; }

    void set_alloc_idx( const size_t& i ){ m_alloc_idx = i; }
    void set_owner( np02_shape *owner ){ m_owner = owner; }
    void set_loc_grid( np02_loc_grid *loc_grid ){m_loc_grid = loc_grid; }
    void set_i( const uint16_t& i){ m_i = i; }
    void set_j( const uint16_t& j){ m_j = j; }
    void set_prev( np02_loc_grid_node *p){ m_prev = p; }
    void set_next( np02_loc_grid_node *n){ m_next = n; }
    void set_s_prev( np02_loc_grid_node *p){ m_s_prev = p; }
    void set_s_next( np02_loc_grid_node *n){ m_s_next = n; }

    const size_t& get_alloc_idx() const{ return m_alloc_idx; }
    np02_shape *get_owner() const { return m_owner; }
    np02_loc_grid *get_loc_grid() const { return m_loc_grid; }
    uint16_t get_i() const { return m_i; }
    uint16_t get_j() const { return m_j; }
    np02_loc_grid_node *get_prev() const { return m_prev; }
    np02_loc_grid_node *get_next() const { return m_next; }
    np02_loc_grid_node *get_s_prev() const { return m_s_prev; }
    np02_loc_grid_node *get_s_next() const { return m_s_next; }
    np02_shp_alloc *get_shp_alloc() const;
    lyr_idx_type get_lyr_idx() const;

    uint64_t hash( const uint64_t& h_in=0 ) const;
    int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    std::ostream& ostream_output(std::ostream& os) const;
    void write_bmp_file(const np02_xy& xy_min,
        const double& pixel_num, const np02_bmp_color& color,
        np02_bmp_file *bmp_file) const;
};

class np02_loc_grid{
public:
    typedef std::vector<np02_loc_grid_node *> loc_grid_node_vec;
    typedef loc_grid_node_vec::iterator loc_grid_node_vec_itr;
    typedef loc_grid_node_vec::const_iterator loc_grid_node_vec_citr;
    typedef std::pair<np02_shape_type, uint32_t> shp_type_idx_pair;
    typedef std::pair<shp_type_idx_pair, np02_shape *> shp_type_idx_shape;
    typedef std::vector<shp_type_idx_shape> shp_type_idx_shape_vec;
    typedef shp_type_idx_shape_vec::const_iterator shp_type_idx_shape_vec_citr;
    typedef shp_type_idx_shape_vec::iterator shp_type_idx_shape_vec_itr;
private:
    np02_shp_alloc *m_shp_alloc;
    size_t m_alloc_idx; /* np02_shp_alloc::alloc_get_loc_grid_by_idx() */
    lyr_idx_type m_lyr_idx;
    np02_loc_grid_dim m_loc_grid_dim;
    double m_extra_search_d;/* extra margin when adding shapes*/
    loc_grid_node_vec m_loc_grid_vec;
    shp_type_idx_shape_vec m_idx_shape_vec;/* temporary vec for shape search */
    np02_uint16_pair_vec m_idx_pair_vec; /* temporary vec */
    double m_small_distance;
public:
    np02_loc_grid();
   ~np02_loc_grid();

    void set_shp_alloc(np02_shp_alloc *shp_alloc){ m_shp_alloc = shp_alloc; }
    void set_alloc_idx( const size_t& i ){ m_alloc_idx = i; }
    void set_lyr_idx( const lyr_idx_type& i ){m_lyr_idx = i; }
    void set_extra_search_d(const double& d){m_extra_search_d = d;}
    void init_loc_grid( const np02_loc_grid_dim& d );
    void insert_shape_in_loc_grid(np02_shape *shape);
    void remove_shape_from_loc_grid(np02_shape *shape);

    np02_shp_alloc *get_shp_alloc() const{ return m_shp_alloc; }
    const size_t& get_alloc_idx() const{ return m_alloc_idx; }
    lyr_idx_type get_lyr_idx() const { return m_lyr_idx; }
    const double& get_extra_search_d() const{ return m_extra_search_d; }
    const np02_loc_grid_dim& get_loc_grid_dim()const{return m_loc_grid_dim;}
    const np02_loc_grid_node *get_loc_grid_head_node(
        const uint16_t i, const uint16_t j){
        const size_t loc_grid_idx = static_cast<size_t>(j) +
           (static_cast<size_t>(i)*static_cast<size_t>(m_loc_grid_dim.get_h()));
        return (loc_grid_idx < m_loc_grid_vec.size()) ?
            m_loc_grid_vec[loc_grid_idx] : NULL; }
    void get_shapes_near_shape(const np02_shape *s,
        np02_shape_vec *shapes) const;
    void get_shapes_in_bb(const np02_xy& xy_min, const np02_xy& xy_max,
        np02_shape_vec *shapes ) const;
    double get_small_distance() const{ return m_small_distance; }

    uint64_t hash( const uint64_t& h_in=0 ) const;
    int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    std::ostream& ostream_output(std::ostream& os) const;
    std::ostream& ostream_output_brief_table(std::ostream& os) const;
    void write_bmp_file(const np02_xy& xy_min,
        const double& pixel_num, const np02_bmp_color& color,
        np02_bmp_file *bmp_file) const;
    void write_bmp_file_grid_lines(const np02_xy& xy_min,
        const double& pixel_num, const np02_bmp_color& color,
        np02_bmp_file *bmp_file) const;
    void write_bmp_file_grid_info(const np02_xy& xy_min,
        const double& pixel_num, const np02_bmp_color& color,
        np02_bmp_file *bmp_file) const;
private:
    np02_loc_grid_node *alloc_loc_grid_node();
    void free_loc_grid_node(np02_loc_grid_node *n);
    void clear_loc_grid();
};

/* Shape Allocator */
class np02_shp_alloc{
private:
    typedef std::vector<np02_circle *> circle_vec;
    typedef std::vector<np02_arc *> arc_vec;
    typedef std::vector<np02_line_seg *> line_seg_vec;
    typedef std::vector<np02_rect *> rect_vec;
    typedef std::vector<np02_polygon *> polygon_vec;
    typedef std::vector<np02_loc_grid_node *> loc_grid_node_vec;
    typedef std::vector<np02_loc_grid *> loc_grid_vec;
private:
    circle_vec m_alloc_circle_vec;
    np02_circle *m_circle_free_chain;
    arc_vec m_alloc_arc_vec;
    np02_arc *m_arc_free_chain;
    line_seg_vec m_alloc_line_seg_vec;
    np02_line_seg *m_line_seg_free_chain;
    rect_vec m_alloc_rect_vec;
    np02_rect *m_rect_free_chain;
    polygon_vec m_alloc_polygon_vec;
    np02_polygon *m_polygon_free_chain;
    loc_grid_node_vec m_alloc_loc_grid_node_vec;
    np02_loc_grid_node *m_loc_grid_node_free_chain;
    loc_grid_vec m_alloc_loc_grid_vec;

public:
    np02_shp_alloc();
    virtual ~np02_shp_alloc();

    void free_shape(np02_shape *shape);

    size_t alloc_get_circle_count() const{ return m_alloc_circle_vec.size(); }
    np02_circle *alloc_get_circle_by_idx(const size_t& i) const{
        return (i<m_alloc_circle_vec.size())?m_alloc_circle_vec.at(i):NULL; }
    np02_circle *alloc_circle();
    void free_circle(np02_circle *circle);

    size_t alloc_get_arc_count() const{ return m_alloc_arc_vec.size(); }
    np02_arc *alloc_get_arc_by_idx(const size_t& i) const{
        return (i<m_alloc_arc_vec.size())?m_alloc_arc_vec.at(i):NULL; }
    np02_arc *alloc_arc();
    void free_arc(np02_arc *arc);

    size_t alloc_get_line_seg_count() const{
        return m_alloc_line_seg_vec.size(); }
    np02_line_seg *alloc_get_line_seg_by_idx(const size_t& i) const{
        return(i<m_alloc_line_seg_vec.size())?m_alloc_line_seg_vec.at(i):NULL;}
    np02_line_seg *alloc_line_seg();
    void free_line_seg(np02_line_seg *line_seg);

    size_t alloc_get_rect_count() const{ return m_alloc_rect_vec.size(); }
    np02_rect *alloc_get_rect_by_idx(const size_t& i) const{
        return (i<m_alloc_rect_vec.size())?m_alloc_rect_vec.at(i):NULL; }
    np02_rect *alloc_rect();
    void free_rect(np02_rect *rect);

    size_t alloc_get_polygon_count() const{return m_alloc_polygon_vec.size();}
    np02_polygon *alloc_get_polygon_by_idx(const size_t& i) const{
        return (i<m_alloc_polygon_vec.size())?m_alloc_polygon_vec.at(i):NULL;}
    np02_polygon *alloc_polygon();
    void free_polygon(np02_polygon *polygon);

    size_t alloc_get_total_shape_count() const;

    size_t alloc_get_loc_grid_node_count() const{
        return m_alloc_loc_grid_node_vec.size(); }
    np02_loc_grid_node *alloc_get_loc_grid_node_by_idx(
        const size_t& i) const{return (i<m_alloc_loc_grid_node_vec.size())?
            m_alloc_loc_grid_node_vec.at(i):NULL; }
    np02_loc_grid_node *alloc_loc_grid_node();
    void free_loc_grid_node(np02_loc_grid_node *loc_grid_node);

    size_t alloc_get_loc_grid_count()const{return m_alloc_loc_grid_vec.size();}
    np02_loc_grid *alloc_get_loc_grid_by_idx(const size_t& i) const{
        return (i<m_alloc_loc_grid_vec.size())?m_alloc_loc_grid_vec.at(i):NULL;}
    np02_loc_grid *alloc_loc_grid();
    void free_loc_grid(np02_loc_grid *loc_grid);

    /* debug */
    uint64_t hash( const uint64_t& h_in=0 ) const;
    virtual int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    virtual std::ostream& ostream_output(std::ostream& os) const;

private:
    void alloc_delete_all();
    int verify_data_alloc_circle( char *err_msg,
        const size_t err_msg_capacity, size_t *err_msg_pos ) const;
    int verify_data_alloc_arc( char *err_msg,
        const size_t err_msg_capacity, size_t *err_msg_pos ) const;
    int verify_data_alloc_line_seg( char *err_msg,
        const size_t err_msg_capacity, size_t *err_msg_pos ) const;
    int verify_data_alloc_rect( char *err_msg,
        const size_t err_msg_capacity, size_t *err_msg_pos ) const;
    int verify_data_alloc_polygon( char *err_msg,
        const size_t err_msg_capacity, size_t *err_msg_pos ) const;
    int verify_data_alloc_loc_grid_node( char *err_msg,
        const size_t err_msg_capacity, size_t *err_msg_pos ) const;
    int verify_data_alloc_loc_grid( char *err_msg,
        const size_t err_msg_capacity, size_t *err_msg_pos ) const;
};

class np02_shape_test{
public:
    static int run_shape_test(const int& shape_test_number, const int& 
        shape_test_iteration_count, const int& shape_test_rand_seed);
private:
    int m_shape_test_number;
    int m_shape_test_iteration_count;
    int m_shape_test_rand_seed;

    /* loop counts */
    uint32_t m_iteration;

    uint32_t m_rand_uint32; /* pseudo-random number*/
public:
    np02_shape_test();
   ~np02_shape_test();
    void set_shape_test_number( const int& shape_test_number ){
        m_shape_test_number = shape_test_number; }
    void set_shape_test_iteration_count(const int& shape_test_iteration_count){
        m_shape_test_iteration_count = shape_test_iteration_count; }
    void set_shape_test_rand_seed( const int& shape_test_rand_seed ){
        m_shape_test_rand_seed = shape_test_rand_seed; }
    int execute();
private:
    void advance_rand();
    int get_rand_int() const { return static_cast<int>(
        static_cast<unsigned int>(
            (m_rand_uint32 >> 16) | (m_rand_uint32 << 16)) & RAND_MAX); }
    double get_rand_dbl(const double& low, const double& high){
        return low + ((high-low)*static_cast<double>(get_rand_int())/
        static_cast<double>(RAND_MAX-1));
        }
    void make_rand_shapes( np02_shp_alloc *shp_alloc, const int& ww,
        const int& hh,  const double& basic_length,
        np02_shape_vec *shapes, std::vector<np02_bmp_color> *colors );
    np02_shape *make_rand_shape( np02_shp_alloc *shp_alloc,
        const np02_xy& shp_ctr, const double& basic_length,
        np02_bmp_color *color );
    void free_shapes( np02_shp_alloc *shp_alloc, np02_shape_vec *shapes );

    int execute_test_1();
    int execute_test_1_0();
    int execute_test_1_1();
    int execute_test_1_2();
    int execute_test_1_i();
    int execute_test_2();
    int execute_test_2_0();
    int execute_test_2_1();
    int execute_test_3();
    int execute_test_4();
    int execute_test_5();
};



class np02_test_main{
public:
    static std::string get_prog_name();
    static std::string get_version_str();
    static int run(int argc, char *argv[]);
private:
    static const char *m_prog_name;
    static const char *m_version_str;
    int m_argc;
    const char **m_argv;

    bool m_test_option;
    int m_test_number;

    bool m_iterations_option;
    int m_iterations;

    bool m_test_iterations_option;
    int m_test_iterations;
    bool m_test_rand_seed_option;
    int m_test_rand_seed;
public:
    np02_test_main(int argc, char *argv[]);
   ~np02_test_main();
    int execute();
private:
    void parse_cmd_line();
    int execute_test_0();
    int execute_test_1();
    int execute_test_2();
    int execute_test_3();
};


void np02_snprintf( char *buf, const size_t buf_capacity, size_t *buf_pos,
               const char *fmt, ... );

} /* namespace np02 */

#endif /* NP02_H */
