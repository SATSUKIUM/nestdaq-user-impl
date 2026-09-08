#ifndef FILTERTIMEFRAMESLICESLICETEMPLATELIB_H
#define FILTERTIMEFRAMESLICESLICETEMPLATELIB_H 1

struct DeleteObject {
  template <typename T>
  void operator()(T *ptr) const
  { delete ptr; }
}; 


#endif // FILTERTIMEFRAMESLICESLICETEMPLATELIB_H
