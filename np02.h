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
class np02_line_seg;
class np02_rect;
class np02_polygon;
class np02_loc_grid_dim;
class np02_loc_grid_node;
class np02_loc_grid;
class np02_shape_test;


class np02_db_object;
class np02_workspace;


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
};

typedef std::vector<np02_xy> np02_xy_vec;
typedef np02_xy_vec::const_iterator np02_xy_citr;
typedef np02_xy_vec::iterator np02_xy_itr;

enum np02_shape_type{
    NP02_SHAPE_TYPE_SHAPE,
    NP02_SHAPE_TYPE_CIRCLE,
    NP02_SHAPE_TYPE_LINE_SEG,
    NP02_SHAPE_TYPE_RECT,
    NP02_SHAPE_TYPE_POLYGON,
    NP02_SHAPE_TYPE_COUNT
};

class np02_shape{
public:
    typedef std::pair<np02_shape_type,uint32_t> type_idx_pair;
    typedef std::pair<type_idx_pair,np02_shape*> type_idx_shp;
    typedef std::vector<type_idx_shp> type_idx_shp_vec;
    typedef type_idx_shp_vec::iterator type_idx_shp_vec_itr;
    typedef type_idx_shp_vec::const_iterator type_idx_shp_vec_citr;
private:
    np02_shape_type m_shape_type; /* type+idx uniquely identifies shape */
    uint32_t m_shape_idx;
    np02_db_object *m_db_obj_owner;
    np02_loc_grid_node *m_head_loc_grid_node;
protected:
    void set_shape_type(const np02_shape_type& shape_type){
        m_shape_type = shape_type; }
public:
    np02_shape();
    virtual ~np02_shape();
    np02_shape_type get_shape_type() const{return m_shape_type;}
    void set_shape_idx(const uint32_t& shape_idx){m_shape_idx=shape_idx;}
    uint32_t get_shape_idx() const{return m_shape_idx;}
    void set_db_obj_owner(np02_db_object *db_obj_owner){
        m_db_obj_owner = db_obj_owner; }
    np02_db_object *get_db_obj_owner() const{return m_db_obj_owner;}
    //np02_layer *get_db_obj_layer() const;
    void set_head_loc_grid_node(np02_loc_grid_node *n){
        m_head_loc_grid_node = n; }
    np02_loc_grid_node *get_head_loc_grid_node() const{
        return m_head_loc_grid_node;}
    np02_loc_grid *get_loc_grid() const;
    np02_workspace *get_workspace() const;
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

    virtual int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    static int verify_distance_from_shape_result(
        const np02_shape *shape_a,
        const np02_shape *shape_b, const np02_xy& xy_near_a,
        const np02_xy& xy_near_b, const double& distance_from,
        char *err_msg, const size_t err_msg_capacity, size_t *err_msg_pos );
    virtual std::ostream& ostream_output(std::ostream& os) const;
    virtual void write_bmp_file(const np02_xy& xy_min,
        const double& pixel_num, const np02_bmp_color& color,
        np02_bmp_file *bmp_file) const;
    virtual void write_dxf_file(const std::string& layer,
        const uint8_t& color, np02_dxf_file *dxf_file) const;
};


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

    virtual int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    virtual std::ostream& ostream_output(std::ostream& os) const;
    virtual void write_bmp_file(const np02_xy& xy_min,
        const double& pixel_num, const np02_bmp_color& color,
        np02_bmp_file *bmp_file) const;
    virtual void write_dxf_file(const std::string& layer,
        const uint8_t& color, np02_dxf_file *dxf_file) const;
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
    virtual double get_distance_from_line_seg_ab(const np02_xy& xy_a,
        const np02_xy& xy_b, np02_xy *near_xy=NULL,
        np02_xy *other_near_xy=NULL) const;
    virtual double get_distance_from_circle(const np02_circle *c,
        np02_xy *near_xy=NULL, np02_xy *circle_near_xy=NULL) const;
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

    int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    std::ostream& ostream_output(std::ostream& os) const;
    void write_bmp_file(const np02_xy& xy_min,
        const double& pixel_num, const np02_bmp_color& color,
        np02_bmp_file *bmp_file) const;
};

/*
  +---------------------------------------------------------------------------------------+
  |np02_shape                                                                         |
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
  |np02_loc_grid                                                                       |
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
    np02_workspace *get_workspace() const;

    int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    std::ostream& ostream_output(std::ostream& os) const;
    void write_bmp_file(const np02_xy& xy_min,
        const double& pixel_num, const np02_bmp_color& color,
        np02_bmp_file *bmp_file) const;
};

class np02_loc_grid{
public:
    typedef std::vector<np02_shape *> shape_vec;
    typedef shape_vec::iterator shape_vec_itr;
    typedef shape_vec::const_iterator shape_vec_citr;
    typedef std::vector<np02_loc_grid_node *> loc_grid_node_vec;
    typedef loc_grid_node_vec::iterator loc_grid_node_vec_itr;
    typedef loc_grid_node_vec::const_iterator loc_grid_node_vec_citr;
    typedef std::pair<np02_shape_type, uint32_t> shp_type_idx_pair;
    typedef std::pair<shp_type_idx_pair, np02_shape *> shp_type_idx_shape;
    typedef std::vector<shp_type_idx_shape> shp_type_idx_shape_vec;
    typedef shp_type_idx_shape_vec::const_iterator shp_type_idx_shape_vec_citr;
    typedef shp_type_idx_shape_vec::iterator shp_type_idx_shape_vec_itr;
private:
    size_t m_alloc_idx;
    //layer *m_owner;
    np02_loc_grid_dim m_loc_grid_dim;
    double m_extra_search_d; /* extra margin when adding shapes to loc*/
    loc_grid_node_vec m_loc_grid_vec;
    shp_type_idx_shape_vec m_idx_shape_vec;/* temporary vec for shape search */
    np02_uint16_pair_vec m_idx_pair_vec; /* temporary vec */
public:
    np02_loc_grid();
   ~np02_loc_grid();

    void set_alloc_idx( const size_t& i ){ m_alloc_idx = i; }
    //void set_owner( layer *owner ){m_owner = owner; }
    void set_extra_search_d(const double& d){m_extra_search_d = d;}
    void init_loc_grid( const np02_loc_grid_dim& d );
    void insert_shape_in_loc_grid(np02_shape *shape);
    void remove_shape_from_loc_grid(np02_shape *shape);

    const size_t& get_alloc_idx() const{ return m_alloc_idx; }
    //layer *get_owner() const { return m_owner; }
    np02_workspace *get_workspace() const;
    const double& get_extra_search_d() const{ return m_extra_search_d; }
    const np02_loc_grid_dim& get_loc_grid_dim()const{return m_loc_grid_dim;}
    const np02_loc_grid_node *get_loc_grid_head_node(
        const uint16_t i, const uint16_t j){
        const size_t loc_grid_idx = static_cast<size_t>(j) +
           (static_cast<size_t>(i)*static_cast<size_t>(m_loc_grid_dim.get_h()));
        return (loc_grid_idx < m_loc_grid_vec.size()) ?
            m_loc_grid_vec[loc_grid_idx] : NULL; }
    void get_shapes_near_shape(const np02_shape *s,
        shape_vec *shapes) const;
    void get_shapes_in_bb(const np02_xy& xy_min, const np02_xy& xy_max,
        shape_vec *shapes ) const;

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


class np02_db_object{
private:
public:
    np02_db_object();
    virtual ~np02_db_object();
    //virtual np02_layer *get_db_obj_layer() const;
    virtual np02_workspace *get_workspace() const=0;
    virtual int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    virtual std::ostream& ostream_output(std::ostream& os) const;
};

class np02_workspace{
private:
    typedef std::vector<np02_circle *> circle_vec;
    typedef std::vector<np02_line_seg *> line_seg_vec;
    typedef std::vector<np02_rect *> rect_vec;
    typedef std::vector<np02_polygon *> polygon_vec;
    typedef std::vector<np02_loc_grid_node *> loc_grid_node_vec;
    typedef std::vector<np02_loc_grid *> loc_grid_vec;
private:
    /* allocator */
    circle_vec m_alloc_circle_vec;
    np02_circle *m_circle_free_chain;
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
    np02_workspace();
    virtual ~np02_workspace();

    /* allocator */
    void free_shape(np02_shape *shape);
    size_t alloc_get_circle_count() const{ return m_alloc_circle_vec.size(); }
    np02_circle *alloc_get_circle_by_idx(const size_t& i) const{
        return (i<m_alloc_circle_vec.size())?m_alloc_circle_vec.at(i):NULL; }
    np02_circle *alloc_circle();
    void free_circle(np02_circle *circle);

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
    virtual int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    virtual std::ostream& ostream_output(std::ostream& os) const;

private:
    void alloc_delete_all();
    int verify_data_alloc_circle( char *err_msg,
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
    void make_rand_shapes( np02_workspace *wksp, const int& ww, 
        const int& hh,  const double& basic_length,
        std::vector<np02_shape *> *shapes,
        std::vector<np02_bmp_color> *colors );
    np02_shape *make_rand_shape( np02_workspace *wksp,
        const np02_xy& shp_ctr, const double& basic_length,
        np02_bmp_color *color );
    void free_shapes( np02_workspace *wksp,
        std::vector<np02_shape *> *shapes );

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
