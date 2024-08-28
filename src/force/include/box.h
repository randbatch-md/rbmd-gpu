#pragma once
#include "types.h"

class Box
{
public:

    enum class BoxType
    {
        ORTHOGONAL,  // 正交
        TRICLINIC    // 三斜
    };

    BoxType _box_type = BoxType::ORTHOGONAL;  // 默认正交盒子

    bool _pbc_x = true;
    bool _pbc_y = true;
    bool _pbc_z = true;


    rbmd::Real* _box_length{0};


};


void MinMirror(Box* box, rbmd::Real& dx, rbmd::Real& dy, rbmd::Real& dz)
{
    //if (box->_box_type == Box::BoxType::ORTHOGONAL)
    //{
    //    if (box->_pbc_x)
    //    {
    //        if (abs(dx > box->_box_length[0] * 0.5))
    //        {
    //            dx -= (dx > 0 ? box->_box_length[0] : -box->_box_length[0]);
    //        }
    //    }
    //    if (box->_pbc_y)
    //    {
    //        if (abs(dy > box->_box_length[1] * 0.5))
    //        {
    //            dy -= (dy > 0 ? box->_box_length[1] : -box->_box_length[1]);
    //        }
    //    }
    //    if (box->_pbc_z)
    //    {
    //        if (abs(dz > box->_box_length[2] * 0.5))
    //        {
    //            dz -= (dz > 0 ? box->_box_length[2] : -box->_box_length[2]);
    //        }
    //    }
    //}
    // TODO else: tri
}
