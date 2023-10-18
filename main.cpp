/* main.cpp  Newport 02 Library

Copyright (c) 2023 Mac Stevens <stevensm@earthlink.net> <www.macstevens.net>

Permission to use, copy, modify, and distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

Reference: https://opensource.org/licenses/ISC
*/
#include "cf01.h"
#include "np02.h"

int main (int argc, char *argv[])
{
//CF01_SET_JRNL_WRITE_MODE_ON();
CF01_SET_JRNL_WRITE_MODE_ON_ERROR();
return np02::np02_test_main::run(argc, argv);
}

