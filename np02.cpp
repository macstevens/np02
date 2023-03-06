/* np02.cpp  Newport 02 Library

Copyright (c) 2023 Mac Stevens <stevensm@earthlink.net> <www.macstevens.net>

Permission to use, copy, modify, and distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

Reference: https://opensource.org/licenses/ISC
*/

#include <algorithm>
#include <iomanip>
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



const char *np02_test_main::m_prog_name = "NP02_TEST";
const char *np02_test_main::m_version_str = "0.0.1" 
     " " __DATE__ " " __TIME__
    ;

std::string np02_test_main::get_prog_name(){ return std::string(m_prog_name); }

std::string np02_test_main::get_version_str(){
char compiler_name[64];
#if defined(__GNUC__)
    sprintf(compiler_name, " GCC %i.%i.%i", __GNUC__, __GNUC_MINOR__,
        __GNUC_PATCHLEVEL__);
#else
    sprintf(compiler_name, "");
#endif
return std::string(m_version_str) 
    //+ std::string(compiler_name)
    ; 
}

int np02_test_main::run(int argc, char *argv[]){
np02_test_main main(argc, argv);
return main.execute();
}

np02_test_main::np02_test_main(int argc, char *argv[]):m_argc(argc),
    m_argv(const_cast<const char **>(&(argv[0]))),
    m_test_option(false),
    m_test_number(0),
    m_test_iterations_option(false),
    m_test_iterations(0),
    m_test_rand_seed_option(false),
    m_test_rand_seed(0)
{
parse_cmd_line();
}

np02_test_main::~np02_test_main()
{}

int np02_test_main::execute(){
int error_code = 0;
if(m_test_option){
    std::cout << m_prog_name << " v" << m_version_str << "\n";
    switch( m_test_number){
        //case 0: error_code = execute_test_0(); break;
        case 1: error_code = execute_test_1(); break;
        case 2: error_code = execute_test_2(); break;
        //case 3: error_code = execute_test_3(); break;
        default:
            std::cout << "attempt to run test " << m_test_number 
                << " fails.  no such test\n";
            break;
        }
    }

if(0 != error_code){
    std::cout << "error code:" << error_code << "\n"; }

return error_code;
}

void np02_test_main::parse_cmd_line(){
/* reset options */
m_test_option = false;
m_test_number = 0;
m_test_iterations_option = false;
m_test_iterations = 0;
m_test_rand_seed_option = false;
m_test_rand_seed = 0;

int i = 1;
bool should_print_help = (m_argc < 2);
while (i < m_argc ){
    if(strncmp(m_argv[i], "--help", strlen("--help"))==0){
        m_test_option=false;
        should_print_help = true;
        ++i;
        }

    else if(strncmp(m_argv[i], "--test=", strlen("--test="))==0){
        m_test_option=true;
        m_test_number = atoi( (m_argv[i]+strlen("--test=")) );
        ++i;
        }
    else if(strncmp(m_argv[i], "-t=", strlen("-t="))==0){
        m_test_option=true;
        m_test_number = atoi( (m_argv[i]+strlen("-t=")) );
        ++i;
        }

    else if(strncmp(m_argv[i], "--iterations=", strlen("--iterations="))==0){
        m_iterations_option=true;
        m_iterations = atoi( (m_argv[i]+strlen("--iterations=")) );
        ++i;
        }
    else if(strncmp(m_argv[i], "--iter=", strlen("--iter="))==0){
        m_iterations_option=true;
        m_iterations = atoi( (m_argv[i]+strlen("--iter=")) );
        ++i;
        }

    else if(strncmp(m_argv[i], "--test-iterations=", strlen("--test-iterations="))==0){
        m_test_iterations_option=true;
        m_test_iterations = atoi( (m_argv[i]+strlen("--test-iterations=")) );
        ++i;
        }
    else if(strncmp(m_argv[i], "--ti=", strlen("--ti="))==0){
        m_test_iterations_option=true;
        m_test_iterations = atoi( (m_argv[i]+strlen("--ti=")) );
        ++i;
        }
    else if(strncmp(m_argv[i], "--test-rand-seed=", strlen("--test-rand-seed="))==0){
        m_test_rand_seed_option=true;
        m_test_rand_seed = atoi( (m_argv[i]+strlen("--test-rand-seed=")) );
        ++i;
        }
    else if(strncmp(m_argv[i], "--seed=", strlen("--seed="))==0){
        m_test_rand_seed_option=true;
        m_test_rand_seed = atoi( (m_argv[i]+strlen("--seed=")) );
        ++i;
        }
    else {
        std::cout << "unknown option: " << m_argv[i] << "\n";
        should_print_help = true;
        ++i;
        }
    }
if(should_print_help){
    std::cout << m_prog_name << " v" << m_version_str << "\n\n";
    std::cout << 
        "Newport 02 Library Test Tool\n\n";

    std::cout << "Regression Tests\n";
    std::cout << "usage: " << m_argv[0] << " --test=<test_number>\n";
    std::cout << "usage: " << m_argv[0] << " --test=0 \\\n"
        "                  --test-rand-seed=<random num gen seed> \\\n"
        "                  --test-iterations=<test iteration count>\n";
    }

/* post-process options*/
if(m_test_option){
    }


}





int np02_test_main::execute_test_1(){
std::cout << "test_1\n";
static const int shape_test_number = 1;
const int shape_test_iteration_count = (m_test_iterations > 0) ?
    m_test_iterations : 1;
const int shape_test_rand_seed = m_test_rand_seed;
const int error_code = np02_shape_test::run_shape_test(shape_test_number, 
        shape_test_iteration_count, shape_test_rand_seed);
return error_code;
}

int np02_test_main::execute_test_2(){
std::cout << "test_2\n";
static const int shape_test_number = 2;
const int shape_test_iteration_count = (m_test_iterations > 0) ?
    m_test_iterations : 1;
const int shape_test_rand_seed = m_test_rand_seed;
const int error_code = np02_shape_test::run_shape_test(shape_test_number, 
        shape_test_iteration_count, shape_test_rand_seed);
return error_code;
}



/* print to buffer
@param buf (output) character buffer
@param buf_capacity (input) size of character buffer
@param buf_pos (input & output)  input: starting position
  output: one character past last non-0x0 character stored
@param fmt C string that contains the text to be written 
*/
void np02_snprintf( char *buf, const size_t buf_capacity, size_t *buf_pos,
               const char *fmt, ... )
{
if((NULL!=buf) &&  (buf_capacity > 0) && (NULL!=buf_pos) &&
    (*buf_pos < buf_capacity) && (NULL != fmt))
    {
    va_list ap;
    size_t n;
    size_t vn;
    va_start(ap, fmt);
    n = buf_capacity - *buf_pos;
    vn = static_cast<size_t>(vsnprintf(buf+*buf_pos, n, fmt, ap));
    *buf_pos = (vn > n) ? buf_capacity : (*buf_pos + vn);
    va_end(ap);
    }
}


} /* namespace np02 */
