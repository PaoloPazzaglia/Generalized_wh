#ifndef __HAS_UNIQUE_ID_HPP__
#define __HAS_UNIQUE_ID_HPP__

template<class X>
class HasUniqueId {
  int id;
  static int counter;
  
  protected:
    void set_id(const int i) {
      id = i;
      counter = (counter > i) ? counter : (i+1);
    }
  
  public:
    HasUniqueId() : id (++counter) {}
    HasUniqueId(const int i) { set_id(i);}
    HasUniqueId(const HasUniqueId& oth) : id(oth.id) { }

    inline int get_id() const { return id; }

};


#endif
