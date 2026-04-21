/*
 *
 *
 */
#ifndef UNPACK_TDC_H
#define UNPACK_TDC_H
#include <iostream>
#include <iomanip>


namespace TDC40 {

static constexpr unsigned int T_TDC = 0xd;
static constexpr unsigned int T_HB = 0xf;
static constexpr unsigned int T_ERROR = 0xe;
static constexpr unsigned int T_SPL_START = 0x1;
static constexpr unsigned int T_SPL_END = 0x4;

struct tdc40 {
	int type;
	int tot;
	int flag;
	int ch;
	int tdc;
	int hartbeat;
};

void Rev5(unsigned char *val, unsigned char *rval);

int Unpack(unsigned char *data, struct tdc40 *tdc);

} // namespace TDC40


namespace TDC64H {

struct tdc64 {
	int type;
	int ch;
	int tot;
	int tdc;
	int tdc4n;
	int flag;
	int spill;
	int hartbeat;
};

static constexpr unsigned int T_TDC        = (0x2c >> 2);
static constexpr unsigned int T_TDC_L      = (0x2c >> 2);
static constexpr unsigned int T_TDC_T      = (0x34 >> 2);
static constexpr unsigned int T_HB         = (0x70 >> 2);
static constexpr unsigned int T_THR1_START = (0x64 >> 2);
static constexpr unsigned int T_THR1_END   = (0x44 >> 2);
static constexpr unsigned int T_THR2_START = (0x68 >> 2);
static constexpr unsigned int T_THR2_END   = (0x48 >> 2);
static constexpr unsigned int T_SPL_START  = (0x60 >> 2);
static constexpr unsigned int T_SPL_END    = (0x50 >> 2);

int Unpack(uint64_t data, struct tdc64 *tdc);

#if 0
int Unpack(uint64_t data, struct tdc64 *tdc)
{
	unsigned char *cdata = reinterpret_cast<unsigned char *>(&data);
	return Unpack(cdata, tdc);
}
#endif

int Unpack(unsigned char *data, struct tdc64 *tdc);

int GetHBFrame(unsigned char *pdata, unsigned char *pend, unsigned char **ppnext);
} //namespace TDC64H

namespace TDC64L {

struct tdc64 {
	int type;
	int ch;
	int tot;
	int tdc;
	int tdc4n;
	int flag;
	int spill;
	int hartbeat;
};

namespace v1 {
static constexpr unsigned int T_TDC       = (0x2c >> 2);
static constexpr unsigned int T_HB        = (0x70 >> 2);
static constexpr unsigned int T_SPL_START = (0x60 >> 2);
static constexpr unsigned int T_SPL_END   = (0x50 >> 2);

int Unpack(uint64_t data, struct tdc64 *tdc);
} //namespace v1

inline namespace v2 {
static constexpr unsigned int T_TDC        = (0x2c >> 2);
static constexpr unsigned int T_TDC_L      = (0x2c >> 2);
static constexpr unsigned int T_TDC_T      = (0x34 >> 2);
static constexpr unsigned int T_THR1_START = (0x64 >> 2);
static constexpr unsigned int T_THR1_END   = (0x44 >> 2);
static constexpr unsigned int T_THR2_START = (0x68 >> 2);
static constexpr unsigned int T_THR2_END   = (0x48 >> 2);
static constexpr unsigned int T_HB         = (0x70 >> 2);
static constexpr unsigned int T_HB1        = (0x70 >> 2);
static constexpr unsigned int T_HB2        = (0x78 >> 2);
static constexpr unsigned int T_SPL_START  = (0x60 >> 2);
static constexpr unsigned int T_SPL_END    = (0x50 >> 2);

int Unpack(uint64_t data, struct tdc64 *tdc);
} //namespace v2

} //namespace TDC64L
 

namespace TDC64H_V3 {

struct tdc64 {
	int type;
	int ch;
	int tot;
	int tdc;
	int tdc4n;
	int flag;
	int toffset;
	int genesize;
	int transize;
	int hartbeat;
};

inline namespace v2 {
static constexpr unsigned int T_TDC        = (0x2c >> 2);
static constexpr unsigned int T_TDC_L      = (0x2c >> 2);
static constexpr unsigned int T_TDC_T      = (0x34 >> 2);
static constexpr unsigned int T_THR1_START = (0x64 >> 2);
static constexpr unsigned int T_THR1_END   = (0x44 >> 2);
static constexpr unsigned int T_THR2_START = (0x68 >> 2);
static constexpr unsigned int T_THR2_END   = (0x48 >> 2);
static constexpr unsigned int T_HB         = (0x70 >> 2);
static constexpr unsigned int T_HB1        = (0x70 >> 2);
static constexpr unsigned int T_HB2        = (0x78 >> 2);
static constexpr unsigned int T_SPL_START  = (0x60 >> 2);
static constexpr unsigned int T_SPL_END    = (0x50 >> 2);

int Unpack(uint64_t data, struct tdc64 *tdc);
} //namespace v2

} //namespace TDC64H_V3
 
namespace TDC64L_V3 {

struct tdc64 {
	int type;
	int ch;
	int tot;
	int tdc;
	int tdc4n;
	int flag;
	int toffset;
	int genesize;
	int transize;
	int hartbeat;
};

inline namespace v2 {
static constexpr unsigned int T_TDC        = (0x2c >> 2);
static constexpr unsigned int T_TDC_L      = (0x2c >> 2);
static constexpr unsigned int T_TDC_T      = (0x34 >> 2);
static constexpr unsigned int T_THR1_START = (0x64 >> 2);
static constexpr unsigned int T_THR1_END   = (0x44 >> 2);
static constexpr unsigned int T_THR2_START = (0x68 >> 2);
static constexpr unsigned int T_THR2_END   = (0x48 >> 2);
static constexpr unsigned int T_HB         = (0x70 >> 2);
static constexpr unsigned int T_HB1        = (0x70 >> 2);
static constexpr unsigned int T_HB2        = (0x78 >> 2);
static constexpr unsigned int T_SPL_START  = (0x60 >> 2);
static constexpr unsigned int T_SPL_END    = (0x50 >> 2);

int Unpack(uint64_t data, struct tdc64 *tdc);
} //namespace v2

} //namespace TDC64L3
 

#ifdef TEST_MAIN_UNPACKTDC
//#if 0

int tdc64h_dump(uint64_t data)
{
	struct TDC64H::tdc64 tdc;
	int type = TDC64H::Unpack(data, &tdc);
	if (type == TDC64H::T_TDC) {
		std::cout << "TDC " << std::dec
			<< " CH:" << tdc.ch 
			<< " TOT:" << tdc.tot 
			<< " TDC:" << tdc.tdc 
			<< " TDC4n:" << tdc.tdc4n 
			<< " : " << std::hex << data
			<< std::endl;
	} else
	if ((type == TDC64H::T_HB)
		|| (tdc.type == TDC64H::T_SPL_START)
		|| (tdc.type == TDC64H::T_SPL_END)) {
		if (type ==TDC64H::T_HB) std::cout << "HB ";
		if (type ==TDC64H::T_SPL_START) std::cout << "S_STA ";
		if (type ==TDC64H::T_SPL_END) std::cout << "S_END ";
		std::cout  << std::hex
			<< " FLAG: " << tdc.flag
			<< " SPILL: " << tdc.spill
			<< " HERTBEAT: " << tdc.hartbeat
			<< " : " << std::hex << data
			<< std::dec << std::endl;
	} else {
		std::cerr << "Invalid data : "
			<< std::hex << data << std::dec << std::endl;
	}

	return type;
}

int main(int argc, char* argv[])
{
	static char cbuf[16];
	uint64_t *pdata = reinterpret_cast<uint64_t *>(cbuf);

	uint64_t *buf = new uint64_t[1024*1024*8];

	int i = 0;
	while (true) {
		std::cin.read(cbuf, 8);
		//std::cin >> *pdata; 
		if (std::cin.eof()) break;
		buf[i++] = *pdata;
		#if 0
		std::cout << "\r "  << i << ": " << *pdata << "  " << std::flush;
		#endif
	
		tdc64h_dump(*pdata);

	}

	unsigned char *pcurr = reinterpret_cast<unsigned char *>(buf);
	unsigned char *pend = reinterpret_cast<unsigned char *>(buf + i);
	unsigned char *pnext = nullptr;
	while (true) {
		int len = TDC64H::GetHBFrame(pcurr, pend, &pnext);
		if (len <= 0) break;

		/*
		std::cout << "#D HB frame size: " << std::dec << len
			<< " curr: " << std::hex
			<< reinterpret_cast<uintptr_t>(pcurr)
			<< " next: "
			<< reinterpret_cast<uintptr_t>(pnext) << std::endl;
		*/

		pdata = reinterpret_cast<uint64_t *>(pcurr);

		for (unsigned int j = 0 ; j < len / sizeof(uint64_t) ; j++) {
			tdc64h_dump(pdata[j]);
		}
		pcurr = pnext;
	}

	return 0;
}
#endif

#if 0
int main(int argc, char* argv[])
{
	static char cbuf[16];
	uint64_t *pdata = reinterpret_cast<uint64_t *>(cbuf);
	struct TDC64H::tdc64 tdc;

	while (true) {
		std::cin.read(cbuf, 8);
		//std::cin >> *pdata; 
		if (std::cin.eof()) break;

		int type = TDC64H::Unpack(*pdata, &tdc);
		if (type == TDC64H::T_TDC) {
			std::cout << "TDC " << std::dec
				<< " CH:" << tdc.ch 
				<< " TOT:" << tdc.tot 
				<< " TDC:" << tdc.tdc 
				<< " : " << std::hex << *pdata
				<< std::endl;
		} else
		if ((type == TDC64H::T_HB)
			|| (tdc.type == TDC64H::T_SPL_START)
			|| (tdc.type == TDC64H::T_SPL_END)) {
			if (type ==TDC64H::T_HB) std::cout << "HB ";
			if (type ==TDC64H::T_SPL_START) std::cout << "S_STA ";
			if (type ==TDC64H::T_SPL_END) std::cout << "S_END ";
			std::cout  << std::hex
				<< " FLAG: " << tdc.flag
				<< " SPILL: " << tdc.spill
				<< " HERTBEAT: " << tdc.hartbeat
				<< " : " << std::hex << *pdata
				<< std::dec << std::endl;
		} else {
			std::cerr << "Invalid data : "
				<< std::hex << *pdata << std::dec << std::endl;
			break;
		}
	}

	return 0;
}
#endif

#if 0
int main(int argc, char* argv[])
{
	static unsigned char cbuf[16];
	//uint64_t *pdata = reinterpret_cast<uint64_t *>(cbuf);
	struct TDC40::tdc40 tdc;

	std::cin.read(reinterpret_cast<char *>(cbuf), 13);
	while (true) {
		std::cin.read(reinterpret_cast<char *>(cbuf), 5);
		//std::cin >> *pdata; 
		if (std::cin.eof()) break;

		int type = TDC40::Unpack(cbuf, &tdc);
		if (type == TDC40::T_TDC) {
			std::cout << "TDC " << std::dec
				<< " CH:" << std::setw(2) << tdc.ch 
				<< " TOT:" << std::setw(3) << tdc.tot 
				<< " TDC:" << std::setw(7) << tdc.tdc 
				<< " : " << std::hex
				<< std::setw(2) << static_cast<unsigned int>(cbuf[4])
				<< std::setw(2) << static_cast<unsigned int>(cbuf[3])
				<< std::setw(2) << static_cast<unsigned int>(cbuf[2])
				<< std::setw(2) << static_cast<unsigned int>(cbuf[1])
				<< std::setw(2) << static_cast<unsigned int>(cbuf[0])
				<< std::endl;
		} else
		if ((type == TDC40::T_HB)
			|| (tdc.type == TDC40::T_ERROR)
			|| (tdc.type == TDC40::T_SPL_START)
			|| (tdc.type == TDC40::T_SPL_END)) {
			if (type ==TDC40::T_HB) std::cout << "HB ";
			if (type ==TDC40::T_ERROR) std::cout << "ERR ";
			if (type ==TDC40::T_SPL_START) std::cout << "S_STA ";
			if (type ==TDC40::T_SPL_END) std::cout << "S_END ";
			std::cout  << std::dec
				//<< " FLAG: " << tdc.flag
				//<< " SPILL: " << tdc.spill
				<< " HERTBEAT: " << tdc.hartbeat
				<< " : " << std::hex
				<< std::setw(2) << static_cast<unsigned int>(cbuf[4])
				<< std::setw(2) << static_cast<unsigned int>(cbuf[3])
				<< std::setw(2) << static_cast<unsigned int>(cbuf[2])
				<< std::setw(2) << static_cast<unsigned int>(cbuf[1])
				<< std::setw(2) << static_cast<unsigned int>(cbuf[0])
				<< std::dec << std::endl;
		} else {
			std::cout << "Invalid data : " << std::hex
				<< std::setw(2) << static_cast<unsigned int>(cbuf[4])
				<< std::setw(2) << static_cast<unsigned int>(cbuf[3])
				<< std::setw(2) << static_cast<unsigned int>(cbuf[2])
				<< std::setw(2) << static_cast<unsigned int>(cbuf[1])
				<< std::setw(2) << static_cast<unsigned int>(cbuf[0])
				<< std::dec << std::endl;
			//break;
		}
	}

	return 0;
}
#endif

#if 0
int main(int argc, char* argv[])
{
	struct TDC40::tdc40 tdc;

	unsigned char rawdata[] = {
		0x00, 0x00, 0x80, 0x40, 0xd0,
		0x11, 0x00, 0x80, 0x40, 0xd0,
		0x22, 0x00, 0x80, 0x40, 0xd0,
		0x33, 0x00, 0x80, 0x40, 0xd0,
		0x44, 0x00, 0x80, 0x40, 0xd0,
		0x02, 0x00, 0x00, 0x00, 0xf0		
	};


	for (int i = 0 ; i < 30 ; i += 5) {
		TDC40::Unpack((rawdata + i), &tdc);
		std::cout << " Head: " << std::hex << tdc.type
			<< " TOT:  " << std::hex << tdc.tot
			<< " TYPE: " << std::hex << tdc.flag
			<< " CH:   " << std::hex << tdc.ch
			<< " TDC:  " << std::hex << tdc.tdc
			<< " HB:  " << std::hex << tdc.hartbeat
			<< std::endl;
	}
	

	return 0;
}
#endif
#endif // UNPACK_TDC_H
