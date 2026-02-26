 #ifndef CLASSMACROS_H
 #define CLASSMACROS_H
 #undef UNCOPYABLE_CLASS
 #if __cplusplus >= 201103L
 #define UNCOPYABLE_CLASS(CLASS_NAME) \
 private:\
 CLASS_NAME(CLASS_NAME const&) = delete;\
 CLASS_NAME& operator=(CLASS_NAME const&) = delete;
 #else
 #define UNCOPYABLE_CLASS(CLASS_NAME) \
 private:\
 CLASS_NAME(CLASS_NAME const&);\
 CLASS_NAME& operator=(CLASS_NAME const&)
 #endif
 #endif
