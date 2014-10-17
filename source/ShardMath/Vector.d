module ShardMath.Vector;
private import std.conv;
private import std.math;
private import std.traits;
import core.stdc.string;
import std.format;
import tested;

// Have to put tests at the top thanks to that @nogc.
@name("Slice Operators")
unittest {
	auto vec = Vector3f(1, 2, 3);
	assert(vec[] == [1, 2, 3]);
	assert(vec[1..2] == [2]);
	vec[0..2] = [2, 3];
	assert(vec == Vector3f(2, 3, 3));
	vec[0..2] = 0;
	assert(vec == Vector3f(0, 0, 3));
}

@name("Unary Operator")
unittest {
	auto vec = Vector3f(2, 3, 4);
	assert(-vec == Vector3f(-2, -3, -4));
	assert(-(-vec) == vec);
}

@name("ToString Operator")
unittest {
	auto vec = Vector3f(1, 2, 3);
	assert(vec.text == "[1, 2, 3]");
	assert(Vector3f(2, -3, 4).text == "[2, -3, 4]");
}

@name("CompareElement Tests")
unittest {
	static if(is(T == float) || is(T == double) || is(T == real)) {
		Vector4f testVector = Vector4f(1.12f, 2.34f, 3.45f, 4.56f);
		assert(testVector.compareElement(0, 1.12f));
		assert(testVector.compareElement(1, 2.34f));
		assert(testVector.compareElement(2, testVector.z));
		assert(testVector.compareElement(3, testVector.w));
		assert(!testVector.compareElement(0, 1));
	}
}

@name("Vector Operators")
unittest {
	auto v1 = Vector3f(1, 2, 3);
	auto v2 = Vector3f(4, 5, 6);
	assert(v1 + v2 == Vector3f(5, 7, 9));
	assert(v2 - v1 == Vector3f(3, 3, 3));
	assert(v2 + v1 + v2 == v2 * 2 + v1);
	assert(v1 + v2 - v2 == v1);
	assert(v1 * v2 == Vector3f(4, 10, 18));
	assert(v1 + Vector3f(0, 0, 0) == v1);
}

@name("Scalar Operators")
unittest {
	auto vec = Vector3f(1, 2, 3);
	assert(vec * 3 == Vector3f(3, 6, 9));
	assert(vec * 3.5f == Vector3f(3.5f, 7f, 10.5f));
	assert(vec / 2 == Vector3f(0.5f, 1f, 1.5f));
	assert(vec + 1 == Vector3f(2, 3, 4));
}

@nogc:

/// A structure representing an N-Dimensional vector of T type.
struct Vector(size_t N, T) if(N >= 2) {	
	
	// To allow outside classes to easier access.
	alias T ElementType;
	alias N NumElements;

	/// If T is a floating point type, this is aliased to T; otherwise, float.
	/// This is used for operations that do not have a meaningful result for integer vectors.
	alias DecimalType = Select!(isFloatingPoint!(T), T, float);

	/// Allows foreach over individual elements of the vector.
	int opApply(int delegate(ref T) dg) {
		int result;
		foreach(ref T element; this.elements)
			if((result = dg(element)) != 0)
				break;
		return result;
	}

	/// Returns this vector as a string in the same way that an array would be.
	void toString(scope void delegate(const(char)[]) sink, FormatSpec!char fmt) const {
		sink("[");
		foreach(i, ref e; elements) {
			sink.formatValue(e, fmt);
			if(i != N - 1)
				sink(", ");
		}
		sink("]");
	}

@safe pure nothrow:
	
	/// Initializes a new instance of the Vector structure.
	/// Params: val = The value to initialize all elements to.
	this(T val) {
		for(size_t i = 0; i < N; i++)
			this.elements[i] = val;
	}
	
	static if(N == 2) {
		/// Initializes a new instance of the Vector structure.			
		this(T x, T y) {
			this.x = x;
			this.y = y;
		}
	} else static if(N == 3) {
		/// Initializes a new instance of the Vector structure.
		this(T x, T y, T z) {
			this.x = x;
			this.y = y;
			this.z = z;
		}
	} else static if(N == 4) {
		/// Initializes a new instance of the Vector structure.
		this(T x, T y, T z, T w) {
			this.x = x; 
			this.y = y;
			this.z = z; 
			this.w = w;
		}
	} else {
		/// Initializes a new instance of the Vector structure.
		this(in T[N] elements...) {
			for(size_t i =0; i < elements.length; i++)	
				this.elements[i] = elements[i];
		}
	}

	/// Basic operator implementations.
	T opIndex(size_t index) const {
		return this.elements[index];
	}

	/// Ditto
	void opIndexAssign(size_t index, T val) {
		this.elements[index] = val;
	}

	/// Ditto
	void opSliceAssign(T val, size_t i1, size_t i2) {
		this.elements[i1..i2] = val;
	}

	/// Ditto
	void opSliceAssign(T[] val, size_t i1, size_t i2) {
		this.elements[i1..i2] = val;
	}
	
	/// Ditto
	T[] opSlice() {
		return elements;
	}
	
	/// Ditto
	T[] opSlice(size_t first, size_t last) {
		// TODO: Return a smaller Vector?
		return elements[first .. last];	
	}

	/// Ditto
	bool opEquals(in Vector!(N, T) other) const {
		for(size_t i = 0; i < N; i++)
			if(!compareElement(i, other.elements[i]))
				return false;
		return true;
	}

	/// Ditto
	bool opEquals(in T[] elements) const {
		return this.elements == elements;
	}

	/// Ditto
	bool opEquals(in ref T[N] elements) const {
		return this.elements == elements;
	}
	
	/// Ditto
	Vector!(N, T) opUnary(string s)() const if(s == "-") {
		Vector!(N, T) result = void;
		for(size_t i = 0; i < N; i++)
			result.elements[i] = -elements[i];
		return result;
	}

	/// Ditto
	void opAssign(in Vector!(N, T) other) @trusted {
		memcpy(elements.ptr, other.elements.ptr, N * T.sizeof);
	}

	/// Binary operators for non-vector scalars.
	Vector opBinary(string Op)(in T scalar) const {
		// TODO: SIMD
		Vector result = void;
		for(size_t i = 0; i < N; i++)
			mixin("result.elements[i] = this.elements[i] " ~ Op ~ " scalar;");
		return result;
	}

	/// Binary operators for multiple vectors.
	Vector opBinary(string Op)(in Vector other) const {
		// TODO: SIMD
		Vector result = void;
		for(size_t i = 0; i < N; i++)
			mixin("result.elements[i] = this.elements[i] " ~ Op ~ " other.elements[i];");
		return result;
	}

	// Unions for how to access the fields.
	static if(N == 2) {
		union {
			struct { T x, y; };
			struct { T width, height; }
			T[2] elements;	
		}
	} else static if(N == 3) {
		union {
			struct { T x, y, z; };			
			T[3] elements;
		}
	} else static if(N == 4) {
		union {
			T[4] elements;
			struct { T x, y, z, w; };			
		}
	} else {		
		T[N] elements;	
	}

private:		
	static if(is(T == float) || is(T == double) || is(T == real)) {
		bool compareElement(size_t index, T value) @safe const pure nothrow {
			return approxEqual(elements[index], value);
		}
	} else {
		bool compareElement(size_t index, T value) @safe const pure nothrow {
			return elements[index] == value;
		}
	}
	
	//static if(((N * T.sizeof) % 16) == 0) {
	//align(16):
	//private enum bool UseAlignedSSE = true;
	//} else
	enum bool UseAlignedSSE = false;
}

@safe pure nothrow:

/// Returns the dot product between the two vectors.
T dot(size_t N, T)(in Vector!(N, T) first, in Vector!(N, T) second) {		
	//TODO: SIMD
	static if(N <= 4) {
		static if(N >= 2) {
			T x = first.x * second.x;
			T y = first.y * second.y;
		}
		static if(N >= 3)
			T z = first.z * second.z;
		static if(N >= 4)
			T w = first.w * second.w;
		static if(N == 2)
			return x + y;
		else static if(N == 3)
			return x + y + z;
		else static if(N == 4)
			return x + y + z + w;
	} else {
		T Sum = 0;
		for(size_t i = 0; i < N; i++)
			Sum += first.elements[i] * second.elements[i];
	}
}

@name("Dot Product Tests")
unittest {
	auto v1 = Vector3f(1, 2, 3);
	auto v2 = Vector3f(4, -5, 6);
	assert(dot(v1, v2) == 12);
}

/// Returns the cross product between the two vectors. Only valid when N is equal to 3.
Vector!(3, T) cross(T)(in Vector!(3, T) first, in Vector!(3, T) second) {
	return Vector!(3, T)(
		first.y * second.z - second.y * first.z, 
		first.z * second.x - second.z * first.x,
		first.x * second.y - second.x * first.y
	);
}

@name("Cross Product Tests")
unittest {
	auto v1 = Vector3f(3, -3, 1);
	auto v2 = Vector3f(4, 9, 2);
	assert(v1.cross(v2) == Vector3f(-15, -2, 39));
}

/// Returns a new Vector with the minimum component of each Vector.
Vector!(N, T) min(size_t N, T)(in Vector!(N, T) first, in Vector!(N, T) second) {
	// TODO: SIMD
	// TODO: Not use reals thanks to std.math.
	static if(N == 2)
		return Vector!(N, T)(cast(T)fmin(first.x, second.x), cast(T)fmin(first.y, second.y));
	else static if(N == 3)
		return Vector!(N, T)(cast(T)fmin(first.x, second.x), cast(T)fmin(first.y, second.y), cast(T)fmin(first.z, second.z));
	else static if(N == 4)
		return Vector!(N, T)(cast(T)fmin(first.x, second.x), cast(T)fmin(first.y, second.y), cast(T)fmin(first.z, second.z), cast(T)fmin(first.w, second.w));
	else {
		T[N] elements;
		for(size_t i = 0; i < N; i++)
			elements[i] = cast(T)fmin(first.elements[i], second.elements[i]);
	}
}

/// Returns a new Vector with the maximum component of each Vector.
Vector!(N, T) max(size_t N, T)(in Vector!(N, T) first, in Vector!(N, T) second) {
	//TODO: Same as for min.
	static if(N == 2)
		return Vector!(N, T)(cast(T)fmax(first.x, second.x),cast(T) fmax(first.y, second.y));
	else static if(N == 3)
		return Vector!(N, T)(cast(T)fmax(first.x, second.x), cast(T)fmax(first.y, second.y), cast(T)fmax(first.z, second.z));
	else static if(N == 4)
		return Vector!(N, T)(cast(T)fmax(first.x, second.x), cast(T)fmax(first.y, second.y), cast(T)fmax(first.z, second.z), cast(T)fmax(first.w, second.w));
	else {
		T[N] elements;
		for(size_t i = 0; i < N; i++)
			elements[i] = cast(T)fmax(first.elements[i], second.elements[i]);
	}
}

@name("Min-Max Tests")
unittest {
	Vector3f v1 = Vector3f(1, 2, 3);
	Vector3f v2 = Vector3f(2, -3, 4);
	assert(v1.min(v2) == Vector3f(1, -3, 3));
	assert(v1.max(v2) == Vector3f(2, 2, 4));
}

/// Returns the distance between the two vectors squared.
T distanceSquared(size_t N, T)(in Vector!(N, T) first, in Vector!(N, T) second) {
	static if(N <= 4) {
		static if(N >= 2) {
			T dX = first.x - second.x;
			T dY = first.y - second.y;
		} 
		static if(N >= 3)
			T dZ = first.z - second.z;
		static if(N >= 4)
			T dW = first.w - second.w;		
	}				
	static if(N == 2)
		return (dX * dX) +(dY * dY);
	else static if(N == 3)
		return (dX * dX) + (dY * dY) + (dZ * dZ);
	else static if(N == 4)
		return (dX * dX) + (dY * dY) + (dZ * dZ) + (dW * dW);
	else {
		T[N] data;
		for(size_t i = 0; i < N; i++) {
			T dI = first.elements[i] - second.elements[i];
			data[i] = dI * dI;				
		}
		return Vector!(N, T)(data);
	}
}

@name("Distance Tests")
unittest {
	auto v1 = Vector2f(3, 2), v2 = Vector2f(9, 7);
	assert(distanceSquared(v1, v2) == 61);
	assert(distanceSquared(Vector2i(3, 2), Vector2i(9, 7)) == 61);
}

/// Linearly interpolates between the first and second vector by the specified amount.
/// Params: Amount = The amount to interpolate by. A value of 0 is entirely the first vector, with 1 being entirely the second vector.	
Vector!(N, T) lerp(size_t N, T)(in Vector!(N, T) first, in Vector!(N, T) second, float Amount) {		
	static if(N <= 4) {
		static if(N >= 2) {
			T X = cast(T)(first.x + ((second.x - first.x) * Amount));
			T Y = cast(T)(first.y + ((second.y - first.y) * Amount));
		}
		static if(N >= 3)
			T Z = cast(T)(first.z + ((second.z - first.z) * Amount));
		static if(N >= 4)
			T W = cast(T)(first.w + ((second.w - first.w) * Amount));
		static if(N == 2)
			return Vector!(N, T)(X, Y);
		else static if(N == 3)
			return Vector!(N, T)(X, Y, Z);
		else
			return Vector!(N, T)(X, Y, Z, W);
	} else {
		T[N] elements;
		for(size_t i = 0; i < N; i++)
			elements[i] = cast(T)(first.elements[i] + ((second.elements[i] - first.elements[i]) * Amount));
		return Vector!(N, T)(elements);
	}
}

@name("Lerp Tests")
unittest {
	auto v1 = Vector3f(1, 2, 3);
	auto v2 = Vector3f(2, 3, 4);
	assert(lerp(v1, v2, 0) == v1);
	assert(lerp(v1, v2, 1) == v2);
	assert(lerp(v1, v2, 0.5f) == Vector3f(1.5f, 2.5f, 3.5f));
}

/// Returns the sum of all of the components in this Vector.
T sum(size_t N, T)(Vector!(N, T) vec) {
	T result = 0;
	foreach(element; vec.elements)
		result += element;
	return result;
}

@name("Sum Tests")
unittest {
	assert(Vector3f(1, 2, 3).sum == 6);
	assert(Vector3f(-1, 0, 1).sum == 0);
}

/// Gets the magnitude, or length, of this Vector.
@property Vector!(N, T).DecimalType magnitude(size_t N, T)(Vector!(N, T) vec) {		
	return cast(Vector!(N, T).DecimalType)sqrt(cast(Vector!(N, T).DecimalType)vec.magnitudeSquared);
}

/// Gets the magnitude, or length, of this Vector without performing a square-root operation on the end result.
@property T magnitudeSquared(size_t N, T)(in Vector!(N, T) vec) {
	T result = 0;
	foreach(ref e; vec.elements)
		result += e * e;
	return result;
}

@name("Magnitude Tests")
unittest {
	auto v = Vector3f(1, 2, 3);
	assert(v.magnitude == sqrt(14f));
	assert(v.magnitudeSquared == 14);
}

/// Returns whether this Vector is normalized (aka, a unit vector).
@property bool isNormalized(size_t N, T)(in Vector!(N, T) vec) {
	return approxEqual(vec.magnitudeSquared - 1, 0);
}

/// Returns a normalized version of this Vector.
Vector!(N, T) normalized(size_t N, T)(in Vector!(N, T) vec) {		
	//TODO: SIMD
	alias DecimalType = vec.DecimalType;
	static if(N <= 4) {
		static if(N >= 2) {
			DecimalType xS = vec.x * vec.x;
			DecimalType yS = vec.y * vec.y;
		} 
		static if(N >= 3)
			DecimalType zS = vec.z * vec.z;
		static if(N >= 4)
			DecimalType wS = vec.w * vec.w;
		static if(N == 2) { 
			auto recip = sqrt(1 / (xS + yS));
			return Vector!(N, T)(cast(T)(vec.x * recip), cast(T)(vec.y * recip));
		} else static if(N == 3) {
			auto recip = sqrt(1 / (xS + yS + zS));
			return Vector!(N, T)(cast(T)(vec.x * recip), cast(T)(vec.y * recip), cast(T)(vec.z * recip));
		} else {
			auto recip = sqrt(1 / (xS + yS + zS + wS));
			return Vector!(N, T)(cast(T)(vec.x * recip), cast(T)(vec.y * recip), cast(T)(vec.z * recip), cast(T)(vec.w * recip));
		}				
	} else {
		DecimalType recip = 0;
		for(size_t i = 0; i < N; i++)
			recip += vec.elements[i] * vec.elements[i];		
		recip = sqrt(1 / recip);
		T[N] data;
		for(size_t i = 0; i < N; i++)
			data = cast(T)(vec.elements[i] * recip);
		return Vector!(N, T)(data);
	}
}

/// Normalizes a Vector2 instance in-place.
void normalizeInline(size_t N, T)(ref Vector!(N, T) vec) {
	alias DecimalType = vec.DecimalType;
	static if(N <= 4) {
		static if(N >= 2) {
			DecimalType xS = vec.x * vec.x;
			DecimalType yS = vec.y * vec.y;
		} 
		static if(N >= 3) {
			DecimalType zS = vec.z * vec.z;
		} 
		static if(N == 4) {
			DecimalType wS = vec.w * vec.w;						
		} 
		static if(N == 2) { 
			auto recip = sqrt(1 / (xS + yS));
			vec.x *= recip; 				
			vec.y *= recip;				
		}  else static if(N == 3) {
			auto recip = sqrt(1 / (xS + yS + zS));
			vec.x *= recip; 				
			vec.y *= recip;				
			vec.z *= recip;		
		} else static if(N == 4) {
			auto recip = sqrt(1 / (xS + yS + zS + wS));
			vec.x *= recip; 				
			vec.y *= recip;				
			vec.z *= recip;		
			vec.w *= recip;				
		}				
	} else {
		DecimalType recip = 0;			
		for(size_t i = 0; i < N; i++)
			recip += vec.elements[i] * vec.elements[i];
		recip = sqrt(1 / recip);
		for(size_t i = 0; i < N; i++)
			vec.elements[i] *= recip;
	}
}

@name("Normalization Tests")
unittest {
	auto v = Vector3f(1, 0, 0);
	assert(v.isNormalized);
	assert(v.normalized == v);
	auto v2 = Vector3f(1, 1, 1);
	assert(!v2.isNormalized);
	auto exNorm = Vector3f(1/sqrt(3f), 1/sqrt(3f), 1/sqrt(3f));
	assert(v2.normalized == exNorm);
	assert(!v2.isNormalized);
	v2.normalizeInline();
	assert(v2 == exNorm);
	assert(v2.isNormalized);
}

///	Determines whether this Vector contains an element with the specified value.
bool contains(size_t N, T)(Vector!(N, T) vec, T val) {	
	for(size_t i = 0; i < N; i++)
		if(vec.compareElement(i, val))
			return true;				
	return false;
}

@name("Contains Tests")
unittest {
	auto v = Vector3f(1, 2, 3);
	assert(v.contains(2));
	assert(!v.contains(4));
}

/// Basic type aliases:
alias Vector!(4, float) Vector4f;
/// Ditto
alias Vector!(4, double) Vector4d;
/// Ditto
alias Vector!(4, int) Vector4i;
/// Ditto
alias Vector!(3, float) Vector3f;
/// Ditto
alias Vector!(3, double) Vector3d;
/// Ditto
alias Vector!(3, int) Vector3i;
/// Ditto
alias Vector!(2, float) Vector2f;
/// Ditto
alias Vector!(2, double) Vector2d;
/// Ditto
alias Vector!(2, int) Vector2i;

//version(Windows)
//	enum bool IsWin32 = true;
//else
enum bool IsWin32 = false;