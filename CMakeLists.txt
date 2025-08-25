cmake_minimum_required(VERSION 4.0)
project(evmlr_shuffle C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(evmlr_shuffle
        evmlr_commit.c
        evmlr_commit.h
        evmlr_enc.c
        evmlr_enc.h
        evmlr_extra.c
        evmlr_extra.h
        evmlr_params.h
        evmlr_utils.c
        evmlr_utils.h)
