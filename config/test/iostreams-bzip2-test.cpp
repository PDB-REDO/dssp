#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

namespace io = boost::iostreams;

int main()
{
	io::filtering_stream<io::input> in;
	in.push(io::bzip2_decompressor());

	return 0;
}

