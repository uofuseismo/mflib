#ifndef MFLIB_PRIVATE_PRIVATE_HPP
#define MFLIB_PRIVATE_PRIVATE_HPP
#ifndef HAVE_STD__MAKE_UNIQUE
#include <memory>
namespace std
{
no
template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
}
#endif
#endif
