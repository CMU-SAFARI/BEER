/**
 * @file gf2.h
 *
 * @brief Representation of GF(2) data types as compatible with Eigen
 *
 * @author Minesh Patel (minesh.patelh@gmail.com)
 */
#ifndef GF2_H
#define GF2_H

#include "Eigen/Eigen"
#include "z3++.h"

/**
 * @brief representation of a GF(2) element
 * 
 * Note: DOES NOT SUPPORT MULTIPLE Z3 CONTEXTS!! IT USES A SINGLE STATIC GLOBAL ONE.
 */
class gf2
{	
private:
	z3::expr e;

public:
	static z3::context *ctx;
	static void init_with_z3_context(z3::context &z3_ctx) { gf2::ctx = &z3_ctx; }

	gf2() : e(*ctx) {}
	gf2(const gf2 &v) : e(v.get_z3_expr()) {}
	gf2(const z3::expr &v) : e(v) {}
	gf2(bool n) : e(ctx->bool_val(n)) {}
	gf2(int n) : e(ctx->bool_val(n)) { assert((n & 1) == n && "creating GF(2) value from non-boolean"); }
	~gf2() {}

	const z3::expr &get_z3_expr(void) const { return e; }
	
	z3::expr &operator=(z3::expr const &n)
	{
		return e = n;
	}

	gf2 &operator=(bool n) { e = ctx->bool_val(n); return *this; }
	gf2 &operator=(int n)  { assert((n & 1) == n && "creating GF(2) value from non-boolean"); e = ctx->bool_val(n); return *this; }

	gf2 operator==(const gf2 &rhs) const { return rhs == this->e; }
	gf2 operator==(const z3::expr &rhs) const { return this->e == rhs; }
	gf2 operator==(const bool rhs) const { return this->e == ctx->bool_val(rhs); }
	
    gf2 operator!=(const gf2 &rhs) const { return rhs != this->e; }
	gf2 operator!=(const z3::expr &rhs) const { return this->e != rhs; }
	gf2 operator!=(const bool rhs) const { return this->e != ctx->bool_val(rhs); }
	
	gf2 operator&&(const gf2 &rhs) const { return rhs && this->e; }
	gf2 operator&&(const z3::expr &rhs) const { return this->e && rhs; }
	gf2 operator&&(const bool rhs) const { return this->e && rhs; }
	
	gf2 operator||(const gf2 &rhs) const { return rhs || this->e; }
	gf2 operator||(const z3::expr &rhs) const { return this->e || rhs; }
	gf2 operator||(const bool rhs) const { return this->e || rhs; }
	
	gf2 operator*(const gf2 &rhs) const { return this->e && rhs.get_z3_expr(); }
	gf2 operator*(const z3::expr &rhs) const { return this->e && rhs; }
	gf2 operator*(const bool rhs) const { return this->e && rhs; }
	
	gf2 operator+(const gf2 &rhs) const { return this->e != rhs.get_z3_expr(); }
	gf2 operator+(const z3::expr &rhs) const { return this->e != rhs; }
	gf2 operator+(const bool rhs) const { return this->e != rhs; }
	
	gf2 &operator+=(const gf2 &rhs) { this->e = this->e != rhs.get_z3_expr(); return *this; }
	gf2 &operator+=(const z3::expr &rhs) { this->e = this->e != rhs; return *this; }
	gf2 &operator+=(const bool rhs) { this->e = this->e != rhs; return *this; }
	
	gf2 operator!(void) const { return !this->e; }

	friend std::ostream &operator<<(std::ostream &os, gf2 const &m);
};

/**
 * @brief Eigen support for processing GF(2) types
 */
namespace Eigen
{
	/**
	 * @brief defining the Eigen traits for the GF(2) type
	 */
	template<> struct NumTraits<gf2>
	    : NumTraits<double> // permits to get the epsilon, dummy_precision, lowest, highest functions
	{
	    typedef gf2 Real; /**< defining GF(2) as a real type */
	    typedef gf2 NonInteger; /**< defining GF(2) as a non-integer type */
	    typedef gf2 Nested; /**< defining GF(2) as a nestable type */
	    enum
	    {
	        IsComplex = 0,
	        IsInteger = 1,
	        IsSigned = 0,
	        RequireInitialization = 1,
	        ReadCost = 1,
	        AddCost = 3,
	        MulCost = 3
	    };
	};

	/**
	 * @brief defining binary operations for GF(2) types on bools
	 */
	template<> 
	struct ScalarBinaryOpTraits< gf2, bool >
	{
		enum { Defined = 1 };
		typedef gf2 ReturnType; /**< defining the return type of the binary operation */
	};

	/**
	 * @brief defining binary operations for bools on GF(2) types
	 * 
	 * @tparam element type
	 */
	template<> 
	struct ScalarBinaryOpTraits< bool, gf2 >
	{
		enum { Defined = 1 };
		typedef gf2 ReturnType; /**< defining the return type of the binary operation */
	};
}


#endif /* GF2_H */
