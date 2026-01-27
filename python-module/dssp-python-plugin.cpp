#include <boost/python.hpp>
#include <Python.h>

#include <dssp.hpp>
#include <memory>

struct statistics_wrapper
{
	int get_residues() { return m_stats.count.residues; }
	int get_chains() { return m_stats.count.chains; }
	int get_SS_bridges() { return m_stats.count.SS_bridges; }
	int get_intra_chain_SS_bridges() { return m_stats.count.intra_chain_SS_bridges; }
	int get_H_bonds() { return m_stats.count.H_bonds; }
	int get_H_bonds_in_antiparallel_bridges() { return m_stats.count.H_bonds_in_antiparallel_bridges; }
	int get_H_bonds_in_parallel_bridges() { return m_stats.count.H_bonds_in_parallel_bridges; }

	dssp::statistics m_stats;
};

enum class ladder_direction_type
{
	parallel,
	antiparallel
};

struct residue_info_wrapper
{
};

class dssp_wrapper
{
  public:
	dssp_wrapper(std::string data, int model_nr = 1, int min_poly_proline_stretch_length = 3, bool calculateSurfaceAccessibility = false)
	{
		struct membuf : public std::streambuf
		{
			membuf(char *text, size_t length)
			{
				this->setg(text, text, text + length);
			}
		} buffer(const_cast<char *>(data.data()), data.length());

		std::istream is(&buffer);

		auto f = cif::pdb::read(is);
		m_dssp = std::make_shared<dssp>(f.front(), model_nr, min_poly_proline_stretch_length, calculateSurfaceAccessibility);
	}

	dssp_wrapper(const dssp_wrapper &) = default;
	dssp_wrapper &operator=(const dssp_wrapper &) = default;

	statistics_wrapper get_statistics()
	{
		statistics_wrapper result{ .m_stats = m_dssp->get_statistics() };
		return result;
	}

	struct residue_info_handle
	{
		residue_info_handle(dssp::residue_info info)
			: m_info(info)
		{
		}

		operator dssp::residue_info &() { return m_info; }

		dssp::residue_info m_info;
	};

	class iterator
	{
	  public:
		using iterator_category = std::forward_iterator_tag;
		using value_type = dssp::residue_info;
		using difference_type = std::ptrdiff_t;
		using pointer = value_type *;
		using reference = residue_info_handle;

		iterator(dssp::residue_info info)
			: m_current(info)
		{
		}

		iterator(const iterator &i) = default;
		iterator &operator=(const iterator &i) = default;

		reference operator*() { return {m_current}; }
		pointer operator->() { return &m_current; }

		iterator &operator++()
		{
			m_current = m_current.next();
			return *this;
		}

		iterator operator++(int)
		{
			auto tmp(*this);
			this->operator++();
			return tmp;
		}

		bool operator==(const iterator &rhs) const
		{
			return m_current == rhs.m_current;
		}

	  private:
		dssp::residue_info m_current;
	};

	auto begin() { return iterator(*m_dssp->begin()); }
	auto end() { return iterator({}); }

	auto get(std::string asym_id, int nr)
	{
		return m_dssp->operator[]({ asym_id, nr });
	}

  private:
	std::shared_ptr<dssp> m_dssp;
};

// shamelessly copied from a stack overflow comment
// https://stackoverflow.com/questions/36485840/wrap-boostoptional-using-boostpython

// Custom exceptions
struct AttributeError : std::exception
{
	[[nodiscard]] const char *what() const noexcept override { return "AttributeError exception"; }
};

struct TypeError : std::exception
{
	[[nodiscard]] const char *what() const noexcept override { return "TypeError exception"; }
};

// Set python exceptions
void translate(const std::exception &e)
{
	if (dynamic_cast<const AttributeError *>(&e))
		PyErr_SetString(PyExc_AttributeError, e.what());
	if (dynamic_cast<const TypeError *>(&e))
		PyErr_SetString(PyExc_TypeError, e.what());
}

struct to_python_optional
{
	static PyObject *convert(const std::optional<float> &obj)
	{
		if (obj)
			return boost::python::incref(boost::python::object(*obj).ptr());
		else
			return boost::python::incref(boost::python::object().ptr());
	}
};

struct to_python_list_of_floats
{
	static PyObject *convert(const std::vector<float> &v)
	{
		boost::python::list list;
		for (auto value : v)
			list.append(value);
		return boost::python::incref(boost::python::object(list).ptr());
	}
};

struct to_python_point
{
	static PyObject *convert(const std::tuple<float, float, float> &v)
	{
		boost::python::dict p;
		p["x"] = std::get<0>(v);
		p["y"] = std::get<1>(v);
		p["z"] = std::get<2>(v);
		return boost::python::incref(boost::python::object(p).ptr());
	}
};

struct to_python_partner
{
	static PyObject *convert(const std::tuple<dssp::residue_info, double> &v)
	{
		boost::python::type_info iv = boost::python::type_id<dssp::residue_info>();
		const boost::python::converter::registration* cv = boost::python::converter::registry::query(iv);
		assert(cv != nullptr);
		if (cv == nullptr)
			throw std::runtime_error("Missing registration");

		if (auto &[ri, e] = v; ri)
		{
			auto c = cv->to_python(&ri);
			return boost::python::incref(boost::python::make_tuple(boost::python::handle<>(c), e).ptr());
		}
		else
			return boost::python::incref(boost::python::make_tuple(boost::python::object(), 0).ptr());
	}
};

struct to_python_bridge_partner
{
	static PyObject *convert(const std::tuple<dssp::residue_info, int, bool> &v)
	{
		if (auto &[ri, nr, parallel] = v; ri)
		{
			boost::python::type_info iv = boost::python::type_id<dssp::residue_info>();
			const boost::python::converter::registration* cv = boost::python::converter::registry::query(iv);
			assert(cv != nullptr);
	
			boost::python::type_info dv = boost::python::type_id<ladder_direction_type>();
			const boost::python::converter::registration* ev = boost::python::converter::registry::query(dv);
			assert(ev != nullptr);

			auto c = cv->to_python(&ri);

			ladder_direction_type direction = parallel ? ladder_direction_type::parallel : ladder_direction_type::antiparallel;
			auto e = ev->to_python(&direction);

			return boost::python::incref(boost::python::make_tuple(
				boost::python::handle<>(c), nr, boost::python::handle<>(e)).ptr());
		}
		else
			return boost::python::incref(boost::python::make_tuple(boost::python::object(), 0, boost::python::object()).ptr());
	}
};

bool test_bond_between_residues(PyObject *a, PyObject *b)
{
	const auto &a_ri = boost::python::extract<dssp::residue_info&>(a);
	const auto &b_ri = boost::python::extract<dssp::residue_info&>(b);

	return test_bond(a_ri, b_ri);
}

BOOST_PYTHON_MODULE(mkdssp)
{
	using namespace boost::python;

	register_exception_translator<AttributeError>(&translate);
	register_exception_translator<TypeError>(&translate);

	to_python_converter<std::optional<float>, to_python_optional>();
	to_python_converter<std::vector<float>, to_python_list_of_floats>();
	to_python_converter<std::tuple<float, float, float>, to_python_point>();
	to_python_converter<std::tuple<dssp::residue_info, double>, to_python_partner>();
	to_python_converter<std::tuple<dssp::residue_info, int, bool>, to_python_bridge_partner>();

	enum_<dssp::structure_type>("structure_type")
		.value("Loop", dssp::structure_type::Loop)
		.value("Alphahelix", dssp::structure_type::Alphahelix)
		.value("Betabridge", dssp::structure_type::Betabridge)
		.value("Strand", dssp::structure_type::Strand)
		.value("Helix_3", dssp::structure_type::Helix_3)
		.value("Helix_5", dssp::structure_type::Helix_5)
		.value("Helix_PPII", dssp::structure_type::Helix_PPII)
		.value("Turn", dssp::structure_type::Turn)
		.value("Bend", dssp::structure_type::Bend);

	enum_<dssp::helix_type>("helix_type")
		.value("_3_10", dssp::helix_type::_3_10)
		.value("alpha", dssp::helix_type::alpha)
		.value("pi", dssp::helix_type::pi)
		.value("pp", dssp::helix_type::pp);

	enum_<dssp::helix_position_type>("helix_position_type")
		.value("NoHelix", dssp::helix_position_type::None)
		.value("Start", dssp::helix_position_type::Start)
		.value("End", dssp::helix_position_type::End)
		.value("StartAndEnd", dssp::helix_position_type::StartAndEnd)
		.value("Middle", dssp::helix_position_type::Middle);

	enum_<ladder_direction_type>("ladder_direction_type")
		.value("parallel", ladder_direction_type::parallel)
		.value("anti_parallel", ladder_direction_type::antiparallel);

	enum_<dssp::chain_break_type>("chain_break_type")
		.value("NoGap", dssp::chain_break_type::None)
		.value("NewChain", dssp::chain_break_type::NewChain)
		.value("Gap", dssp::chain_break_type::Gap);

	class_<statistics_wrapper>("statistic_counts")
		.add_property("residues", &statistics_wrapper::get_residues)
		.add_property("chains", &statistics_wrapper::get_chains)
		.add_property("SS_bridges", &statistics_wrapper::get_SS_bridges)
		.add_property("intra_chain_SS_bridges", &statistics_wrapper::get_intra_chain_SS_bridges)
		.add_property("H_bonds", &statistics_wrapper::get_H_bonds)
		.add_property("H_bonds_in_antiparallel_bridges", &statistics_wrapper::get_H_bonds_in_antiparallel_bridges)
		.add_property("H_bonds_in_parallel_bridges", &statistics_wrapper::get_H_bonds_in_parallel_bridges);

	class_<dssp::residue_info>("residue_info")
		.add_property("asym_id", &dssp::residue_info::asym_id)
		.add_property("seq_id", &dssp::residue_info::seq_id)
		.add_property("alt_id", &dssp::residue_info::alt_id)
		.add_property("compound_id", &dssp::residue_info::compound_id)
		.add_property("compound_letter", &dssp::residue_info::compound_letter)
		.add_property("auth_asym_id", &dssp::residue_info::auth_asym_id)
		.add_property("auth_seq_id", &dssp::residue_info::auth_seq_id)
		.add_property("pdb_strand_id", &dssp::residue_info::pdb_strand_id)
		.add_property("pdb_seq_num", &dssp::residue_info::pdb_seq_num)
		.add_property("pdb_ins_code", &dssp::residue_info::pdb_ins_code)
		.add_property("alpha", &dssp::residue_info::alpha)
		.add_property("kappa", &dssp::residue_info::kappa)
		.add_property("phi", &dssp::residue_info::phi)
		.add_property("psi", &dssp::residue_info::psi)
		.add_property("tco", &dssp::residue_info::tco)
		.add_property("omega", &dssp::residue_info::omega)
		.add_property("is_pre_pro", &dssp::residue_info::is_pre_pro)
		.add_property("is_cis", &dssp::residue_info::is_cis)
		.add_property("chiral_volume", &dssp::residue_info::chiral_volume)
		.add_property("chi", &dssp::residue_info::chis)
		.add_property("ca_location", &dssp::residue_info::ca_location)
		.add_property("chain_break", &dssp::residue_info::chain_break)
		.add_property("nr", &dssp::residue_info::nr)
		.add_property("type", &dssp::residue_info::type)
		.add_property("ssBridgeNr", &dssp::residue_info::ssBridgeNr)
		.def("helix", &dssp::residue_info::helix, args("helix_type"), "Return the position of this residue in a helix with type helix_type")
		.add_property("is_alpha_helix_end_before_start", &dssp::residue_info::is_alpha_helix_end_before_start)
		.add_property("bend", &dssp::residue_info::bend)
		.add_property("accessibility", &dssp::residue_info::accessibility)
		.def("bridge_partner", &dssp::residue_info::bridge_partner, args("indexnr"), "Return a tuple containing the residue, number and direction for the bridge partner with index indexnr")
		.add_property("sheet", &dssp::residue_info::sheet)
		.add_property("strand", &dssp::residue_info::strand)
		.def("acceptor", &dssp::residue_info::acceptor, args("indexnr"), "Return a tuple containing the residue and bond energy for the acceptor with index indexnr")
		.def("donor", &dssp::residue_info::donor, args("indexnr"), "Return a tuple containing the residue and bond energy for the donor with index indexnr");

	class_<dssp_wrapper, boost::noncopyable>("dssp", init<std::string, optional<int, int, bool>>())
		.add_property("statistics", &dssp_wrapper::get_statistics)
		.def("__iter__", iterator<dssp_wrapper>())
		.def("get", &dssp_wrapper::get, args("asym_id", "seq_id"), "Return the residue info object for the residue with specified asym_id and seq_id");

	def("TestBond", test_bond_between_residues, args("a", "b"), "Returns true if residues a and b are bonded according to DSSP");
}