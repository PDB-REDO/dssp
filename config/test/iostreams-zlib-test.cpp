#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace io = boost::iostreams;

int main()
{
	io::filtering_stream<io::input> in;
	in.push(io::gzip_decompressor());

	return 0;
}

