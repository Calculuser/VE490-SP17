#ifndef _MDARRAY_H_
#define _MDARRAY_H_

#include <vector>
#include <memory.h>

template <class T>
class CMDArray
{
public:
  explicit CMDArray(int s1, int s2 = -1, int s3 = -1, int s4 = -1, int s5 = -1, int s6 = -1, int s7 = -1, int s8 = -1, int s9 = -1, int s10 = -1, int s11 = -1, int s12 = -1) : data_(NULL)
  {
    data_ = NULL;
    Reset(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12);
  }

  CMDArray()
  {
    data_ = NULL;
  }

  bool Reset(int s1, int s2 = -1, int s3 = -1, int s4 = -1, int s5 = -1, int s6 = -1, int s7 = -1, int s8 = -1, int s9 = -1, int s10 = -1, int s11 = -1, int s12 = -1)
  {
    if (data_ != NULL)
    {
      delete[] data_;
      data_ = NULL;
    }

    if (s1 <= 0)
    {
      data_ = NULL;
      return false;
    }

    memset(sizes_, 0, sizeof(sizes_));
    memset(bases_, 0, sizeof(bases_));

    sizes_[0] = s1;
    bases_[0] = 1;
    dimension_count_ = 1;

    if (s2 > 0)
    {
      sizes_[1] = s2;
      bases_[1] = bases_[0] *  sizes_[0];
      dimension_count_ = 2;
    }

    if (s3 > 0)
    {
      sizes_[2] = s3;
      bases_[2] = bases_[1] *  sizes_[1];
      dimension_count_ = 3;
    }

    if (s4 > 0)
    {
      sizes_[3] = s4;
      bases_[3] = bases_[2] * sizes_[2];
      dimension_count_ = 4;
    }

    if (s5 > 0)
    {
      sizes_[4] = s5;
      bases_[4] = bases_[3] * sizes_[3];
      dimension_count_ = 5;
    }

    if (s6 > 0)
    {
      sizes_[5] = s6;
      bases_[5] = bases_[4] * sizes_[4];
      dimension_count_ = 6;
    }

    if (s7 > 0)
    {
      sizes_[6] = s7;
      bases_[6] = bases_[5] * sizes_[5];
      dimension_count_ = 7;
    }

    if (s8 > 0)
    {
      sizes_[7] = s8;
      bases_[7] = bases_[6] * sizes_[6];
      dimension_count_ = 8;
    }

    if (s9 > 0)
    {
      sizes_[8] = s9;
      bases_[8] = bases_[7] * sizes_[7];
      dimension_count_ = 9;
    }

    if (s10 > 0)
    {
      sizes_[9] = s10;
      bases_[9] = bases_[8] * sizes_[8];
      dimension_count_ = 10;
    }

    if (s11 > 0)
    {
      sizes_[10] = s11;
      bases_[10] = bases_[9] * sizes_[9];
      dimension_count_ = 11;
    }

    if (s12 > 0)
    {
      sizes_[11] = s12;
      bases_[11] = bases_[10] * sizes_[10];
      dimension_count_ = 12;
    }

    int total_size = bases_[dimension_count_ - 1] * sizes_[dimension_count_ - 1];
    data_ = new T[total_size];
      memset(data_, 0, sizeof(T) * total_size);
    return true;
  }

  virtual ~CMDArray()
  {
    delete[] data_;
    data_ = NULL;
  }

  T& operator()(int i1)
  {
    return data_[bases_[0] * i1];
  }

  T& operator()(int i1, int i2)
  {
    return data_[bases_[0] * i1 + bases_[1] * i2];
  }

  T& operator()(int i1, int i2, int i3)
  {
    return data_[bases_[0] * i1 + bases_[1] * i2 + bases_[2] * i3];
  }

  T& operator()(int i1, int i2, int i3, int i4)
  {
    return data_[bases_[0] * i1 + bases_[1] * i2 + bases_[2] * i3 + bases_[3] * i4];
  }

  T& operator()(int i1, int i2, int i3, int i4, int i5)
  {
    return data_[bases_[0] * i1 + bases_[1] * i2 + bases_[2] * i3 + bases_[3] * i4 + bases_[4] * i5];
  }

  T& operator()(int i1, int i2, int i3, int i4, int i5, int i6)
  {
    return data_[bases_[0] * i1 + bases_[1] * i2 + bases_[2] * i3 + bases_[3] * i4 + bases_[4] * i5 + bases_[5] * i6];
  }

  T& operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7)
  {
    return data_[bases_[0] * i1 + bases_[1] * i2 + bases_[2] * i3 + bases_[3] * i4 + bases_[4] * i5 + bases_[5] * i6 + bases_[6] * i7];
  }

  T& operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8)
  {
    return data_[bases_[0] * i1 + bases_[1] * i2 + bases_[2] * i3 + bases_[3] * i4 + bases_[4] * i5 + bases_[5] * i6 + bases_[6] * i7 + bases_[7] * i8];
  }

  T& operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9)
  {
    return data_[bases_[0] * i1 + bases_[1] * i2 + bases_[2] * i3 + bases_[3] * i4 + bases_[4] * i5 + bases_[5] * i6 + bases_[6] * i7 + bases_[7] * i8, bases_[8] * i9];
  }

  T& operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, int i10)
  {
    return data_[bases_[0] * i1 + bases_[1] * i2 + bases_[2] * i3 + bases_[3] * i4 + bases_[4] * i5 + bases_[5] * i6 + bases_[6] * i7 + bases_[7] * i8, bases_[8] * i9 + bases_[9] * i10];
  }

  T& operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, int i10, int i11)
  {
    return data_[bases_[0] * i1 + bases_[1] * i2 + bases_[2] * i3 + bases_[3] * i4 + bases_[4] * i5 + bases_[5] * i6 + bases_[6] * i7 + bases_[7] * i8, bases_[8] * i9 + bases_[9] * i10 + bases_[10] * i11];
  }

  T& operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, int i10, int i11, int i12)
  {
    return data_[bases_[0] * i1 + bases_[1] * i2 + bases_[2] * i3 + bases_[3] * i4 + bases_[4] * i5 + bases_[5] * i6 + bases_[6] * i7 + bases_[7] * i8, bases_[8] * i9 + bases_[9] * i10 + bases_[10] * i11 + bases_[11] * i12];
  }

  //T& operator()(int i1, int i2 = -1, int i3 = -1, int i4 = -1, int i5 = -1, int i6 = -1, int i7 = -1, int i8 = -1, int i9 = -1, int i10 = -1, int i11 = -1, int i12 = -1)
  //{
  //  size_t offset = 0;

  //  switch (dimension_count_)
  //  {
  //  case 12:
  //    offset += bases_[11] * i12;
  //  case 11:
  //    offset += bases_[10] * i11;
  //  case 10:
  //    offset += bases_[9] * i10;
  //  case 9:
  //    offset += bases_[8] * i9;
  //  case 8:
  //    offset += bases_[7] * i8;
  //  case 7:
  //    offset += bases_[6] * i7;
  //  case 6:
  //    offset += bases_[5] * i6;
  //  case 5:
  //    offset += bases_[4] * i5;
  //  case 4:
  //    offset += bases_[3] * i4;
  //  case 3:
  //    offset += bases_[2] * i3;
  //  case 2:
  //    offset += bases_[1] * i2;
  //  case 1:
  //    offset += bases_[0] * i1;
  //  }

  //  return data_[offset];
  //}

  void set_zero()
  {
    if (data_ != NULL)
    {
      int total_size = bases_[dimension_count_ - 1] * sizes_[dimension_count_ - 1];
      memset(data_, 0, sizeof(T) * total_size);
    }
  }

private:
  T* data_;
  int sizes_[12];
  int bases_[12];
  size_t dimension_count_;
};


#endif

