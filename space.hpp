#ifndef _SPACE_HPP_
#define _SPACE_HPP_
#include <array>

// class Matrix
// {
  
// };

template<unsigned... index_set_shape>
class FiniteSet
{
  
};

// class Topology : public Set
// {
  
// };

template<int N>
class EuclideanVectorSpaceElement
{
  std::array<double,N> element;
  double inner_product_();
};


class Vector : public EuclideanVectorSpaceElement
{
  // Vector add_(){};
  Vector  operator+(const Vector& added_v);
  // Vector subtract_(){};
  Vector operator-(const Vector& subtract_v);
  // Vector scalar_multiply_(){};
  friend Vector operator*(const double& scalar,const Vector& v_self);
};


class EuclideanSpacePoint
{
  std::array<double,3> point;
  EuclideanSpacePoint operator+(const Vector& added_v);
  EuclideanSpacePoint operator-(const Vector& added_v);
};


template <class T> class FieldAtPoint
{
  EuclideanSpacePoint point;
  T field_value;
};



class CoordinateSystem : public FiniteSet
{
  std::array<int,index_set_shape> element;
  Vector
};


#endif
