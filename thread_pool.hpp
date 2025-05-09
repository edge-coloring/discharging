#include <boost/asio.hpp>
#include <boost/thread.hpp>
#include <boost/bind/bind.hpp>
#include <boost/shared_ptr.hpp>

class ThreadPool {
    boost::asio::io_service& io_service_;
    boost::shared_ptr<boost::asio::io_service::work> work_;
    boost::thread_group group_;
public:
    ThreadPool(boost::asio::io_service& io_service, size_t size): io_service_(io_service) {
        // work を設定しておくことで実行しているスレッドの数が 0 になってもスレッドプールが動き続けるようにする。
        work_.reset(new boost::asio::io_service::work(io_service_));
        // size の分だけスレッドを作り、io_service_.run しておく。
        for (size_t i = 0; i < size; ++i) {
            group_.create_thread(boost::bind(&boost::asio::io_service::run, &io_service_));
        }
    }

    ~ThreadPool() {
        work_.reset();
        group_.join_all();
    }

    template <class F>
    void post(F f) {
        io_service_.post(f);
        return;
    }
};
