/*
 * Copyright (c) 2018 Subbaraya Sundeep <sundeep.babi@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
*/
#include "emu.h"

#define LOG(fmt, ...)				\
do {	\
	if (log)	\
		fprintf(stderr, "%x: %x "fmt" %d\n", cpu.pc, cpu.ins,	\
			## __VA_ARGS__, cpu.ins_count);	\
} while(0)

#define TEST(x, y) ((x) & (1 << (y)))
#define ARRAY_SIZE(x) (sizeof(x)/sizeof x[0])

#define SB(b, n) ((b) |= (1 << (n)))
#define CB(b, n) ((b) &= ~(1 << (n)))

#define SBCF(b)	((b) |= (1 << (0)))
#define CBCF(b)	((b) &= ~(1 << (0)))

#define SBZF(b)	((b) |= (1 << (1)))
#define CBZF(b)	((b) &= ~(1 << (1)))

#define SBNF(b)	((b) |= (1 << (2)))
#define CBNF(b)	((b) &= ~(1 << (2)))

#define SBVF(b)	((b) |= (1 << (3)))
#define CBVF(b)	((b) &= ~(1 << (3)))

static int log;

static inline uint16_t extract16(uint16_t value, int start, int length)
{
	return (value >> start) & (0xFFFF >> (16 - length));
}

char ins_names[][10] = {
	"LDX", "STX", "LDX+", "STX+", "LD-X", "ST-X", "PUSH", "POP", "LDS", "STS",
	"LDY+", "STY+", "LD-Y", "ST-Y", "LDZ+", "STZ+", "LD-Z", "ST-Z", "LPMZ",
	"ELPMZ", "LPMZ+", "ELPMZ+", "XCH", "LAS", "LAC", "LAT", "NEG", "SWAP", "INC",
	"ASR", "LSR", "ROR", "SEX", "CLX", "RET", "RETI", "SLEEP", "BREAK", "WDR",
	"LPM", "ELPM", "SPM", "SPM1", "INDJMP", "DEC", "DES", "JMP", "CALL", "ADIW",
	"SBIW", "CBI", "SBI", "SBIC", "SBIS", "MUL", "OUT", "IN", "RJMP", "RCALL",
	"LDI", "BRANCH", "BLD", "BST", "SBRS", "SBRC", "LDDY", "LDDZ", "STDY", "STDZ",
	"NOP", "MOVW", "MULS", "MULSU", "FMUL", "FMULS", "CPSE", "CPC", "SBC",
	"ADC", "AND", "EOR", "OR", "MOV", "CPI", "SBCI", "ORI", "ANDI"
};

char ioreg_names[][15] = {
	"TWBR", "TWSR", "TWAR", "TWDR", "ADCL", "ADCH", "ADCSRA", "ADMUX", "ACSR",
	"UBRRL", "UCSRB", "UCSRA", "UDR", "SPCR", "SPSR", "SPDR", "PIND", "DDRD",
	"PORTD", "PINC", "DDRC", "PORTC", "PINB", "DDRB", "PORTB", "Reserved1",
	"Reserved2", "Reserved3", "EECR", "EEDR", "EEARL", "EEARH", "UCSRC",
	"WDTCR", "ASSR", "OCR2", "TCNT2", "TCCR2", "ICR1L", "ICR1H", "OCR1BL",
	"OCR1BH", "OCR1AL", "OCR1AH", "TCNT1L", "TCNT1H", "TCCR1B", "TCCR1A",
	"SFIOR", "OSCCAL", "TCNT0", "TCCR0", "MCUCSR", "MCUCR", "TWCR", "SPMCR",
	"TIFR", "TIMSK", "GIFR", "GICR", "Reserved4", "SPL", "SPH", "SREG"
};

static uint16_t bp_table[20];

enum instructions {
	LDX = 0,
	STX,
	LDXi,
	STXi,
	LDdX,
	STdX,
	PUSH,
	POP,
	LDS,
	STS,
	LDYi,
	STYi,
	LDdY,
	STdY,
	LDZi,
	STZi,
	LDdZ,
	STdZ,
	LPMZ,
	ELPMZ,
	LPMZi,
	ELPMZi,
	XCH,
	LAS,
	LAC,
	LAT,
	NEG,
	SWAP,
	INC,
	ASR,
	LSR,
	ROR,
	RET,
	RETI,
	SLEEP,
	BREAK,
	WDR,
	LPM,
	ELPM,
	SPM,
	SPMZi,
	SECL,
	INDJMP,
	DEC,
	DES,
	JMP,
	CALL,
	ADIW,
	SBIW,
	CBI,
	SBI,
	SBIC,
	SBIS,
	MUL,
	OUT,
	IN,
	RJMP,
	RCALL,
	LDI,
	BRANCH,
	BLD,
	BST,
	SBRS,
	SBRC,
	LDDY,
	LDDZ,
	STDY,
	STDZ,
	NOP,
	MOVW,
	MULS,
	MULSU,
	FMUL,
	FMULS,
	CPSE,
	CPC,
	CP,
	SBC,
	SUB,
	ADC,
	ADD,
	AND,
	EOR,
	OR,
	MOV,
	CPI,
	SBCI,
	SUBI,
	ORI,
	ANDI,
	BKPT,
};

enum ioregs {
	TWBR = 0x0,
	TWSR,
	TWAR,
	TWDR,
	ADCL,
	ADCH,
	ADCSRA,
	ADMUX,
	ACSR,
	UBRRL,
	UCSRB,
	UCSRA,
	UDR,
	SPCR,
	SPSR,
	SPDR,
	PIND,
	DDRD,
	PORTD,
	PINC,
	DDRC,
	PORTC,
	PINB,
	DDRB,
	PORTB,
	Reserved1,
	Reserved2,
	Reserved3,
	EECR,
	EEDR,
	EEARL,
	EEARH,
	UCSRC,
	WDTCR,
	ASSR,
	OCR2,
	TCNT2,
	TCCR2,
	ICR1L,
	ICR1H,
	OCR1BL,
	OCR1BH,
	OCR1AL,
	OCR1AH,
	TCNT1L,
	TCNT1H,
	TCCR1B,
	TCCR1A,
	SFIOR,
	OSCCAL,
	TCNT0,
	TCCR0,
	MCUCSR,
	MCUCR,
	TWCR,
	SPMCR,
	TIFR,
	TIMSK,
	GIFR,
	GICR,
	Reserved4,
	SPL,
	SPH,
	SREG,
};

#define SP ((cpu.sram.regs.ioreg[SPH] << 8) | cpu.sram.regs.ioreg[SPL])

struct registers {
	uint8_t reg[32]; /* registers r0-r31 */
	uint8_t ioreg[64]; /* io registers */
} __attribute__((packed));

union data_mem {
	struct registers regs;  
	uint8_t ram[1120];
} __attribute__((packed));

struct avr {
	uint16_t pc:12;
	uint16_t ins;
	uint8_t rom[8192];
	union data_mem sram;
	int ins_idx;
	uint8_t rd;
	uint8_t rr;
	uint8_t io;
	uint16_t simm;
	uint16_t imm;
	uint8_t bn:3;
	uint8_t bv:1;
	uint8_t bitnum:3;
	/* 0: free_run , 1: under_monitor */
	bool mode;
	uint32_t ins_count;
} cpu;

void io_wr(uint8_t io, uint8_t val)
{
	switch(io) {
	case UDR:
		if (cpu.sram.regs.ioreg[UCSRB] & (1 << 3))
			write(1, &val, 1);
		break;
	default:
		break;
	}
}

uint8_t io_rd(uint8_t io)
{
	uint8_t val = 0;

	switch(io) {
	case UDR:
		if (cpu.sram.regs.ioreg[UCSRB] & (1 << 4))
			fscanf(stdin, "%c", &val);
		break;
	default:
		val = cpu.sram.regs.ioreg[io];		
		break;
	}

	return val;
}

static void store_ram(uint16_t addr, uint8_t byte)
{
	if (addr > sizeof(union data_mem)) {
		printf("Invalid addr:0x%x\n", addr);
		exit(1);
	}

	if ((addr >= 0x20) && (addr <= 0x5F)) {
		io_wr(addr - 0x20, byte);
	}
	cpu.sram.ram[addr] = byte;
}

static uint8_t load_ram(uint16_t addr)
{
	uint8_t byte;

	if (addr > sizeof(union data_mem)) {
		printf("Invalid addr:0x%x\n", addr);
		exit(1);
	}

	if ((addr >= 0x20) && (addr <= 0x5F)) {
		return io_rd(addr - 0x20);
	}
	return cpu.sram.ram[addr]; 
}

static inline void push(uint8_t val)
{
	uint16_t tmp;

	tmp = SP;
	store_ram(tmp, val);
	tmp--;
	store_ram(SPL + 0x20, tmp & 0xFF);
	store_ram(SPH + 0x20, (tmp & 0xFF00) >> 8);
}

static inline uint8_t pop(void)
{
	uint16_t tmp;

	tmp = SP;
	tmp ++;
	store_ram(SPL + 0x20, tmp & 0xFF);
	store_ram(SPH + 0x20, (tmp & 0xFF00) >> 8);

	return load_ram(tmp);
}

static inline void XYZ(char ch, char inc, uint16_t val, bool flags)
{
	uint16_t tmp;
	uint8_t r1, r2;

	switch (ch) {
	case 'W':
		r1 = 24;
		r2 = 25;
		break;
	case 'X':
		r1 = 26;
		r2 = 27;
		break;
	case 'Y':
		r1 = 28;
		r2 = 29;
		break;
	case 'Z':
		r1 = 30;
		r2 = 31;
		break;
	default:
		return;
		break;
	}

	tmp = (cpu.sram.regs.reg[r2] << 8) | cpu.sram.regs.reg[r1];
	if (inc == '+')
		tmp += val;
	else
		tmp -= val;
	cpu.sram.regs.reg[r1] = tmp;
	cpu.sram.regs.reg[r2] = tmp >> 8;
}

static inline uint16_t W(void)
{
	return (cpu.sram.regs.reg[25] << 8) | cpu.sram.regs.reg[24];
}
static inline uint16_t X(void)
{
	return (cpu.sram.regs.reg[27] << 8) | cpu.sram.regs.reg[26];
}
static inline uint16_t Y(void)
{
	return (cpu.sram.regs.reg[29] << 8) | cpu.sram.regs.reg[28];
}
static inline uint16_t Z(void)
{
	return (cpu.sram.regs.reg[31] << 8) | cpu.sram.regs.reg[30];
}

uint32_t loadbin()
{
	FILE *fp;
	uint32_t len;

	fp = fopen("./uart.bin", "rb");
	if (!fp) {
		printf("file error");
		return 0;
	}

	fseek(fp, 0, SEEK_END);
	len = ftell(fp);
	rewind(fp);

	fread(cpu.rom, len, 2, fp);
	fclose(fp);

	return len; 
}

static uint16_t fetch(uint16_t pc)
{
	return *(uint16_t *)(cpu.rom + cpu.pc);
}

#if 0
static void update_flags(void)
{
	cpu.sram.regs.ioreg[SREG] = (cpu.I << 7) | (cpu.T << 6) | (cpu.H << 5) |
		(cpu.S << 4) | (cpu.V << 3) | (cpu.N << 2) | (cpu.Z << 1) | cpu.C;
}
#endif

static uint8_t add(uint8_t r1, uint8_t r2, bool C)
{
	uint16_t tmp;
	uint8_t ret;
	uint16_t tmp1 = r1;
	uint16_t tmp2 = r2;
	uint8_t *sreg;
	uint8_t carry = 0;

	sreg = &cpu.sram.regs.ioreg[SREG];
	if (C && TEST(*sreg, 0))
		carry = 1;

	tmp = ret = tmp1 + tmp2 + carry;

	*sreg &= 0xE0; /* Clear S V N Z C */

	/* TODO: Half Carry
	if ((r1 & 0x08) && (r2 & 0x08))
		cpu.H = 1;
	else
		cpu.H = 0;
	*/

	/* I T H S V N Z C */
	if (tmp > 0xff)
		SB(*sreg, 0);

	if (!ret)
		SB(*sreg, 1);

	if (tmp & 0x80)
		SB(*sreg, 2);
	/*
		ADDITION SIGN BITS
		num1sign num2sign sumsign
		---------------------------
		0 0 0
*OVER*	0 0 1 (adding two positives should be positive)
		0 1 0
		0 1 1
		1 0 0
		1 0 1
*OVER*	1 1 0 (adding two negatives should be negative)
		1 1 1
	*/

	if (!(r1 & 0x80) && !(r2 & 0x80) && (tmp & 0x80))
		SB(*sreg, 3);

	if ((r1 & 0x80) && (r2 & 0x80) && !(tmp & 0x80))
		SB(*sreg, 3);

	tmp = *sreg;
	tmp = extract16(tmp, 2, 1) ^ extract16(tmp, 3, 1);
	if (tmp)
		SB(*sreg, 4);
	else
		CB(*sreg, 4);

	return ret;
}

static uint8_t sub(uint8_t r1, uint8_t r2, bool C)
{
	uint8_t ret;
	uint8_t *sreg;
	uint8_t carry = 0;

	sreg = &cpu.sram.regs.ioreg[SREG];
	if (C && TEST(*sreg, 0))
		carry = 1;

	ret = r1 - r2 - carry;

	*sreg &= 0xE0; /* Clear S V N Z C */

	/* TODO: Half carry
	if ((r1 & 0x08) && (r2 & 0x08))
		cpu.H = 1;
	else
		cpu.H = 0;
	*/
	/* I T H S V N Z C */
	if ((r2 + carry) > r1)
		SB(*sreg, 0);

	if (!ret)
		SB(*sreg, 1);

	if (ret & 0x80)
		SB(*sreg, 2);
	/*
		SUBTRACTION SIGN BITS
		num1sign num2sign sumsign
		---------------------------
		0 0 0
		0 0 1
		0 1 0
*OVER*	0 1 1 (subtracting a negative is the same as adding a positive)
*OVER*	1 0 0 (subtracting a positive is the same as adding a negative)
		1 0 1
		1 1 0
		1 1 1
	*/
	if (!(r1 & 0x80) && (r2 & 0x80) && (ret & 0x80))
		SB(*sreg, 3);

	if ((r1 & 0x80) && !(r2 & 0x80) && !(ret & 0x80))
		SB(*sreg, 3);

	if ((*sreg >> 2) ^ (*sreg >> 3))
		SB(*sreg, 4);
	else
		CB(*sreg, 4);

	return ret;
}

struct {
	const char *set;
	const char *clear;
} br[] = {
	{"BRCS", "BRCC"}, /* C */
	{"BREQ", "BRNE"}, /* Z */
	{"BRMI", "BRPL"}, /* N */
	{"BRVS", "BRVC"}, /* V */
	{"BRLT", "BRGE"}, /* S */
	{"BRHS", "BRHC"}, /* H */
	{"BRTS", "BRTC"}, /* T */
	{"BRIE", "BRID"}  /* I */
};

static void handle_br(uint8_t bn, uint8_t bv, uint8_t simm, bool exe)
{
	char ch = ' ';
	uint8_t sreg = cpu.sram.regs.ioreg[SREG];

	/* bit value is inverted in instruction */
	bv = ~bv & 0x1;
	if (simm & 0x40) { /* + or - */
		simm  = (~simm + 1) & 0x7F;
		simm = 2 * simm;
		ch = '-';
	} else {
		simm = 2 * simm;
		ch = '+';
	}
	if (!exe) {
		LOG("%s .%c%d", bv ? br[bn].set : br[bn].clear, ch, simm);
		return;
	}
	if (bv && TEST(sreg, bn))
		goto condition_met;
	else if (!bv && !TEST(sreg, bn))
		goto condition_met;
	else
		return; /* do nothing */

condition_met:
	if (ch == '+')
		cpu.pc += simm;
	else
		cpu.pc -= simm;
}

struct bktbl {
	uint16_t addr;
	uint16_t ins;
	uint16_t immed;
}; 

static struct bktbl bk_table[20];
static int num_bkpts;

static void insert_bkpt(uint16_t addr, uint16_t ins)
{
	bk_table[num_bkpts].addr = addr;	
	bk_table[num_bkpts].ins = ins;
	num_bkpts++;
}


uint16_t get_ins_bptbl(uint16_t addr)
{
	int i;

	for (i = 0; i < num_bkpts; i++)
	{
		if (bk_table[i].addr == addr)
			return bk_table[i].ins;
	}

}

static int execute()
{
	uint16_t tmp;

	switch(cpu.ins_idx) {
	case LDX:
		cpu.sram.regs.reg[cpu.rd] = load_ram(X());
		break;
	case STX:
		store_ram(X(), cpu.sram.regs.reg[cpu.rd]);
		break;
	case LDXi:
		cpu.sram.regs.reg[cpu.rd] = load_ram(X());
		XYZ('X', '+', 1, false);
		break;
	case STXi:
		store_ram(X(), cpu.sram.regs.reg[cpu.rd]);
		XYZ('X', '+', 1, false);
		break;
	case LDdX:
		XYZ('X', '-', 1, false);
		cpu.sram.regs.reg[cpu.rd] = load_ram(X());
		break;
	case STdX:
		XYZ('X', '-', 1, false);
		store_ram(X(), cpu.sram.regs.reg[cpu.rd]);
		break;
	case PUSH:
		push(cpu.sram.regs.reg[cpu.rr]);	
		break;
	case POP:
		cpu.sram.regs.reg[cpu.rr] = pop();	
		break;
	case LDS:
		cpu.sram.regs.reg[cpu.rd] = load_ram(cpu.imm);	
		break;
	case STS:
		store_ram(cpu.imm, cpu.sram.regs.reg[cpu.rd]);
		break;
	case LDYi:
		cpu.sram.regs.reg[cpu.rd] = load_ram(Y());
		XYZ('Y', '+', 1, false);
		break;
	case STYi:
		store_ram(Y(), cpu.sram.regs.reg[cpu.rd]);
		XYZ('Y', '+', 1, false);
		break;
	case LDdY:
		XYZ('Y', '-', 1, false);
		cpu.sram.regs.reg[cpu.rd] = load_ram(Y());
		break;
	case STdY:
		XYZ('Y', '-', 1, false);
		store_ram(Y(), cpu.sram.regs.reg[cpu.rd]);
		break;
	case LDZi:
		cpu.sram.regs.reg[cpu.rd] = load_ram(Z());
		XYZ('Z', '+', 1, false);
		break;
	case STZi:
		store_ram(Z(), cpu.sram.regs.reg[cpu.rd]);
		XYZ('Z', '+', 1, false);
		break;
	case LDdZ:
		XYZ('Z', '-', 1, false);
		cpu.sram.regs.reg[cpu.rd] = load_ram(Z());
		break;
	case STdZ:
		XYZ('Z', '-', 1, false);
		store_ram(Z(), cpu.sram.regs.reg[cpu.rd]);
		break;
	case LPMZ:
		cpu.sram.regs.reg[cpu.rd] = cpu.rom[Z()];
		break;
	case ELPMZ:
		break;
	case LPMZi:
		cpu.sram.regs.reg[cpu.rd] = cpu.rom[Z()];
		XYZ('Z', '+', 1, false);
		break;
	case ELPMZi:
		break;
	case XCH:
		break;
	case LAS:
		break;
	case LAC:
		break;
	case LAT:
		break;
	case NEG:
		break;
	case SWAP:
		break;
	case INC:
		cpu.sram.regs.reg[cpu.rd] =
		add(cpu.sram.regs.reg[cpu.rd], 1, false);
		break;
	case ASR:
		break;
	case LSR:
		break;
	case ROR:
		break;
	case RET:
		cpu.pc = pop() & 0xFF;
		cpu.pc |= ((pop() & 0xFF) << 8);
		break;
	case RETI:
		cpu.pc = pop() & 0xFF;
		cpu.pc |= ((pop() & 0xFF) << 8);
		SB(cpu.sram.regs.ioreg[SREG], 7);
		break;
	case SLEEP:
		break;
	case BREAK:
		break;
	case WDR:
		break;
	case LPM:
		/* rd is 0 */
		cpu.sram.regs.reg[0] = cpu.rom[Z()];
		break;
	case ELPM:
		break;
	case SPM:
		break;
	case SPMZi:
		break;
	case SECL:
		cpu.bv ? CB(cpu.sram.regs.ioreg[SREG], cpu.bn) : 
			SB(cpu.sram.regs.ioreg[SREG], cpu.bn);
		break;
	case INDJMP:
		cpu.pc = Z();
		break;
	case DEC:
		cpu.sram.regs.reg[cpu.rd] =
		sub(cpu.sram.regs.reg[cpu.rd], 1, false);
		break;
	case DES:
		break;
	case JMP:
		break;
	case CALL:
		break;
	case ADIW:
		if (cpu.rd == 24)
			XYZ('W', '+', cpu.imm, true);
		else if (cpu.rd == 26)
			XYZ('X', '+', cpu.imm, true);
		else if (cpu.rd == 28)
			XYZ('Y', '+', cpu.imm, true);
		else if (cpu.rd == 30)
			XYZ('Z', '+', cpu.imm, true);
		break;
	case SBIW:
		if (cpu.rd == 24)
			XYZ('W', '-', cpu.imm, true);
		else if (cpu.rd == 26)
			XYZ('X', '-', cpu.imm, true);
		else if (cpu.rd == 28)
			XYZ('Y', '-', cpu.imm, true);
		else if (cpu.rd == 30)
			XYZ('Z', '-', cpu.imm, true);
		break;
	case CBI:
		cpu.sram.regs.ioreg[cpu.io] &= (1 << cpu.bitnum);
		break;
	case SBI:
		cpu.sram.regs.ioreg[cpu.io] |= (1 << cpu.bitnum);
		break;
	case SBIC:
		if ((cpu.sram.regs.ioreg[cpu.io] & (1 << cpu.bitnum)) == 0)
			cpu.pc += 2;
		break;
	case SBIS:
		if (cpu.sram.regs.ioreg[cpu.io] & (1 << cpu.bitnum))
			cpu.pc += 2;
		break;
	case MUL:
		break;
	case OUT:
		store_ram(cpu.io + 0x20, cpu.sram.regs.reg[cpu.rd]);
		break;
	case IN:
		cpu.sram.regs.reg[cpu.rd] = load_ram(cpu.io + 0x20);
		break;
	case RJMP:
		if (cpu.simm & 0x0800) { /* + or - */
			tmp = (~cpu.simm + 1) & 0x0FFF;
			tmp = 2 * tmp;
			cpu.pc -= tmp;
		} else {
			tmp = 2 * cpu.simm;
			cpu.pc += tmp;
		}
		break;
	case RCALL:
		push((cpu.pc & 0xFF00) >> 8);
		push(cpu.pc & 0xFF);
		if (cpu.simm & 0x0800) { /* + or - */
			tmp = (~cpu.simm + 1) & 0x0FFF;
			tmp = 2 * tmp;
			cpu.pc -= tmp;
		} else {
			tmp = 2 * cpu.simm;
			cpu.pc += tmp;
		}
		break;
	case LDI:
		cpu.sram.regs.reg[cpu.rd] = cpu.imm;
		break;
	case BRANCH:
		handle_br(cpu.bn, cpu.bv, cpu.simm, true);
		break;
	case BLD:
		break;
	case BST:
		break;
	case SBRS:
		break;
	case SBRC:
		break;
	case LDDY:
		tmp = Y() + cpu.imm;
		cpu.sram.regs.reg[cpu.rd] = load_ram(tmp);
		break;
	case LDDZ:
		tmp = Z() + cpu.imm;
		cpu.sram.regs.reg[cpu.rd] = load_ram(tmp);
		break;
	case STDY:
		tmp = Y() + cpu.imm;
		store_ram(tmp, cpu.sram.regs.reg[cpu.rd]);
		break;
	case STDZ:
		tmp = Z() + cpu.imm;
		store_ram(tmp, cpu.sram.regs.reg[cpu.rd]);
		break;
	case NOP:
		break;
	case MOVW:
		cpu.sram.regs.reg[cpu.rd + 1] = cpu.sram.regs.reg[cpu.rr + 1];
		cpu.sram.regs.reg[cpu.rd] = cpu.sram.regs.reg[cpu.rr];
		break;
	case MULS:
		break;
	case MULSU:
		break;
	case FMUL:
		break;
	case FMULS:
		break;
	case CPSE:
		break;
	case CPC:
		sub(cpu.sram.regs.reg[cpu.rd], cpu.sram.regs.reg[cpu.rr], true);
		break;
	case CP:
		sub(cpu.sram.regs.reg[cpu.rd], cpu.sram.regs.reg[cpu.rr], false);
		break;
	case SBC:
		cpu.sram.regs.reg[cpu.rd] =
		sub(cpu.sram.regs.reg[cpu.rd], cpu.sram.regs.reg[cpu.rr], true);
		break;
	case SUB:
		cpu.sram.regs.reg[cpu.rd] =
		sub(cpu.sram.regs.reg[cpu.rd], cpu.sram.regs.reg[cpu.rr], false);
		break;
	case ADC:
		cpu.sram.regs.reg[cpu.rd] =
		add(cpu.sram.regs.reg[cpu.rd], cpu.sram.regs.reg[cpu.rr], true);
		break;
	case ADD:
		cpu.sram.regs.reg[cpu.rd] =
		add(cpu.sram.regs.reg[cpu.rd], cpu.sram.regs.reg[cpu.rr], false);
		break;
	case AND:
		cpu.sram.regs.reg[cpu.rd] = cpu.sram.regs.reg[cpu.rd] &
									cpu.sram.regs.reg[cpu.rr];
		break;
	case EOR:
		cpu.sram.regs.reg[cpu.rd] = cpu.sram.regs.reg[cpu.rd] ^
									cpu.sram.regs.reg[cpu.rr];
		break;
	case OR:
		cpu.sram.regs.reg[cpu.rd] = cpu.sram.regs.reg[cpu.rd] |
									cpu.sram.regs.reg[cpu.rr];
		break;
	case MOV:
		cpu.sram.regs.reg[cpu.rd] = cpu.sram.regs.reg[cpu.rr];
		break;
	case CPI:
		sub(cpu.sram.regs.reg[cpu.rd], cpu.imm, false);
		break;
	case SBCI:
		cpu.sram.regs.reg[cpu.rd] =
		sub(cpu.sram.regs.reg[cpu.rd], cpu.imm, true);
		break;
	case SUBI:
		cpu.sram.regs.reg[cpu.rd] =
		sub(cpu.sram.regs.reg[cpu.rd], cpu.imm, false);
		break;
	case ORI:
		cpu.sram.regs.reg[cpu.rd] |= cpu.imm;
		if (cpu.sram.regs.reg[cpu.rd])
			CBZF(cpu.sram.regs.ioreg[SREG]);
		else
			SBZF(cpu.sram.regs.ioreg[SREG]);

		if (cpu.sram.regs.reg[cpu.rd] & 0x80)
			SBNF(cpu.sram.regs.ioreg[SREG]);
		else
			CBNF(cpu.sram.regs.ioreg[SREG]);

		CBVF(cpu.sram.regs.ioreg[SREG]);

		if ((cpu.sram.regs.ioreg[SREG] >> 2) ^
				(cpu.sram.regs.ioreg[SREG] >> 3))
			SB(cpu.sram.regs.ioreg[SREG], 4);
		else
			CB(cpu.sram.regs.ioreg[SREG], 4);
		break;
	case ANDI:
		cpu.sram.regs.reg[cpu.rd] &= cpu.imm;
		if (cpu.sram.regs.reg[cpu.rd])
			CBZF(cpu.sram.regs.ioreg[SREG]);
		else
			SBZF(cpu.sram.regs.ioreg[SREG]);

		if (cpu.sram.regs.reg[cpu.rd] & 0x80)
			SBNF(cpu.sram.regs.ioreg[SREG]);
		else
			CBNF(cpu.sram.regs.ioreg[SREG]);

		CBVF(cpu.sram.regs.ioreg[SREG]);

		if ((cpu.sram.regs.ioreg[SREG] >> 2) ^
				(cpu.sram.regs.ioreg[SREG] >> 3))
			SB(cpu.sram.regs.ioreg[SREG], 4);
		else
			CB(cpu.sram.regs.ioreg[SREG], 4);
		break;
	case BKPT:
		cpu.pc -= 2;
		*(uint16_t *)(cpu.rom + cpu.pc) = get_ins_bptbl(cpu.pc);
		return BKPT;
		break;
	default:
		LOG("UNDEFINED");
		break;
	}
	return 0;
}

char flags[] = {'C','Z','N','V','S','H','T','I'};

static void decode_others(uint16_t ins)
{
	uint16_t simm;
	uint8_t rd;
	uint8_t rr;
	uint8_t io;
	uint16_t tmp;

	switch (ins & 0xFC00) {
	case 0x9000:
		tmp = extract16(ins, 0, 4);
		cpu.rd = extract16(ins, 4, 5);
		switch (tmp) {
		case 0x0:
			cpu.imm = fetch(cpu.pc + 2); /* not instruction 16 bit SRAM address */
			if (TEST(ins, 9)) {
				cpu.ins_idx = STS;
				LOG("STS");
			} else {
				cpu.ins_idx = LDS;
				LOG("LDS");
			}
			break;
		case 0x1:
			if (TEST(ins, 9)) {
				cpu.ins_idx = STZi;
				LOG("STZ+ r%d", cpu.rd);
			} else {
				cpu.ins_idx = LDZi;
				LOG("LDZ+ r%d", cpu.rd);
			}
			break;
		case 0x2:
			if (TEST(ins, 9)) {
				cpu.ins_idx = STdZ;
				LOG("ST-Z r%d", cpu.rd);
			} else {
				cpu.ins_idx = LDdZ;
				LOG("LD-Z r%d", cpu.rd);
			}
			break;
		case 0x4:
			if (TEST(ins, 9)) {
				LOG("XCH");
			} else {
				cpu.ins_idx = LPMZ;
				LOG("LPM r%d, Z", cpu.rd);
			}
			break;
		case 0x5:
			if (TEST(ins, 9)) {
				LOG("LAS");
			} else {
				cpu.ins_idx = LPMZi;
				LOG("LPM r%d, Z+", cpu.rd);
			}
			break;
		case 0x6:
			if (TEST(ins, 9)) {
				LOG("LAC");
			} else {
				LOG("ELPM r%d, Z", cpu.rd);
			}
			break;
		case 0x7:
			if (TEST(ins, 9)) {
				LOG("LAT");
			} else {
				LOG("ELPM r%d, Z+", cpu.rd);
			}
			break;
		case 0x9:
			if (TEST(ins, 9)) {
				cpu.ins_idx = STYi;
				LOG("STY+ r%d", cpu.rd);
			} else {
				cpu.ins_idx = LDYi;
				LOG("LDY+ r%d", cpu.rd);
			}
			break;
		case 0xA:
			if (TEST(ins, 9)) {
				cpu.ins_idx = STdY;
				LOG("ST-Y r%d", cpu.rd);
			} else {
				cpu.ins_idx = LDdY;
				LOG("LD-Y r%d", cpu.rd);
			}
			break;
		case 0xC:
			if (TEST(ins, 9)) {
				cpu.ins_idx = STX;
				LOG("STX r%d", cpu.rd);
			} else {
				cpu.ins_idx = LDX;
				LOG("LDX r%d", cpu.rd);
			}
			break;
		case 0xD:
			if (TEST(ins, 9)) {
				cpu.ins_idx = STXi;
				LOG("STX+ r%d", cpu.rd);
			} else {
				cpu.ins_idx = LDXi;
				LOG("LDX+ r%d", cpu.rd);
			}
			break;
		case 0xE:
			if (TEST(ins, 9)) {
				cpu.ins_idx = STdX;
				LOG("ST-X r%d", cpu.rd);
			} else {
				cpu.ins_idx = LDdX;
				LOG("LD-X r%d", cpu.rd);
			}
			break;
		case 0xF:
			rr = extract16(ins, 4, 5);
			cpu.rr = rr;
			if (TEST(ins, 9)) {
				LOG("PUSH r%d", rr);
				cpu.ins_idx = PUSH;
			} else {
				LOG("POP r%d", rr);
				cpu.ins_idx = POP;
			}
			break;
		default:
			break;
		};
		break;
	case 0x9400:
		rr = extract16(ins, 4, 2);
		cpu.imm = ((ins & 0x00C0) >> 2)|(ins & 0xF);
		switch (rr) {
		case 0:
			cpu.rd = 24;
			break;
		case 1:
			cpu.rd = 26;
			break;
		case 2:
			cpu.rd = 28;
			break;
		case 3:
			cpu.rd = 30;
			break;
		default:
			cpu.rd = 24;
			break;
		}
		if (TEST(ins, 9)) {
			if (TEST(ins, 8)) {
				LOG("SBIW r%d, 0x%x", cpu.rd, cpu.imm );
				cpu.ins_idx = SBIW;
			} else {
				LOG("ADIW r%d, 0x%x", cpu.rd, cpu.imm );
				cpu.ins_idx = ADIW;
			}
			break;
		}
		switch (ins & 0xF) {
		case 0: 
			LOG("COM");
			break;
		case 1:
			LOG("NEG");
			break;
		case 2:
			LOG("SWAP");
			break;
		case 3:
			cpu.rd = extract16(ins, 4, 5);
			LOG("INC r%d", cpu.rd);
			cpu.ins_idx = INC;
			break;
		case 5:
			LOG("ASR");
			break;
		case 6:
			LOG("LSR");
			break;
		case 7:
			LOG("ROR");
			break;
		case 8:
			if (TEST(ins, 8)) {
				switch (ins & 0xF0) {
				case 0x00:
					LOG("RET");
					cpu.ins_idx = RET;
					break;
				case 0x10:
					LOG("RETI");
					cpu.ins_idx = RETI;
					break;
				case 0x80:
					LOG("SLEEP");
					break;
				case 0x90:
					LOG("BREAK");
					break;
				case 0xA0:
					LOG("WDR");
					break;
				case 0xC0:
					LOG("LPM");
					break;
				case 0xD0:
					LOG("ELPM");
					break;
				case 0xE0:
					LOG("SPM");
					break;
				case 0xF0:
					LOG("SPMZ+");
					break;
				default:
					break;
				}
			} else {
				cpu.bn = extract16(ins, 4, 3);
				cpu.bv = extract16(ins, 7, 1);
				cpu.ins_idx = SECL;
				LOG("%s%c", cpu.bv ? "CL" : "SE", flags[cpu.bn]);
			}
			break;
		case 9:
			LOG("Indirect jump");
			cpu.ins_idx = INDJMP;
			break;
		case 10:
			cpu.rd = extract16(ins, 4, 5);
			LOG("DEC r%d", cpu.rd);
			cpu.ins_idx = DEC;
			break;
		case 11:
			LOG("DES");
			break;
		case 12: 
		case 13:
			LOG("JMP");
			break;
		case 14:
		case 15:
			LOG("CALL");
			break;
		default:
			break;
		}
		break;
	case 0x9800:
		cpu.bitnum = extract16(ins, 0, 3);
		cpu.io = extract16(ins, 3, 5);
		if (TEST(ins, 8)) {
			LOG("SBIC/SBIS");
			if (TEST(ins, 9)) {
				LOG("SBIC 0x%x, %d\n", cpu.io, cpu.bitnum);
				cpu.ins_idx = SBIC;
			} else {
				LOG("SBIS 0x%x, %d\n", cpu.io, cpu.bitnum);
				cpu.ins_idx = SBIS;
			}
		} else {
			if (TEST(ins, 9)) {
				LOG("SBI 0x%x, %d\n", cpu.io, cpu.bitnum);
				cpu.ins_idx = SBI;
			} else {
				LOG("CBI 0x%x, %d\n", cpu.io, cpu.bitnum);
				cpu.ins_idx = CBI;
			}
		}	
		break;
	case 0x9C00:
		LOG("MUL");
		break;
	default:
		{
			if ((ins & 0xF000) == 0xB000) {
				rd = extract16(ins, 4, 5);
				io = extract16(ins, 0, 4);
				io |= (extract16(ins, 9, 2) << 4);
				cpu.rd = rd;
				cpu.io = io;
				if (TEST(ins, 11)) {
					LOG("OUT 0x%x, r%d", io, rd);
					cpu.ins_idx = OUT;
				}
				else {
					LOG("IN r%d, 0x%x", rd, io);
					cpu.ins_idx = IN;
				}
			}
			if ((ins & 0xF000) == 0xC000) {
				simm = ins & 0x0FFF;
				cpu.simm = simm;
				cpu.ins_idx = RJMP;
				if (cpu.simm & 0x0800) { /* + or - */
					tmp = (~cpu.simm + 1) & 0x0FFF;
					tmp = 2 * tmp;
					LOG("RJMP .-%d", tmp);
				} else {
					tmp = 2 * cpu.simm;
					LOG("RJMP .+%d", tmp);
				}
			}
			if ((ins & 0xF000) == 0xD000) {
				simm = ins & 0x0FFF;
				cpu.simm = simm;
				cpu.ins_idx = RCALL;
				if (cpu.simm & 0x0800) { /* + or - */
					tmp = (~cpu.simm + 1) & 0x0FFF;
					tmp = 2 * tmp;
					LOG("RCALL .-%d", tmp);
				} else {
					tmp = 2 * cpu.simm;
					LOG("RCALL .+%d", tmp);
				}
			}
			if ((ins & 0xF000) == 0xE000) {
				rd = extract16(ins, 4, 4) + 16;
				tmp = (extract16(ins, 8, 4) << 4) | extract16(ins, 0, 4);
				cpu.imm = tmp;
				cpu.rd = rd;
				cpu.ins_idx = LDI;
				LOG("LDI r%d,0x%x", rd, tmp);
			}
			if ((ins & 0xF000) == 0xF000) {
				cpu.bn = extract16(ins, 0, 3);
				cpu.rd = extract16(ins, 4, 5);
				switch (ins & 0x0C00) {
				case 0x0:
				case 0x0400: /* fall through */
					cpu.bv = extract16(ins, 10, 1);
					cpu.simm = extract16(ins, 3, 7);
					handle_br(cpu.bn, cpu.bv, cpu.simm, false);
					cpu.ins_idx = BRANCH;
					break;
				case 0x0800:
					if (TEST(ins, 9)) {
						LOG("BST r%d, r%d", cpu.rd, cpu.bn);
						cpu.ins_idx = BST;
					} else {
						LOG("BLD r%d, r%d", cpu.rd, cpu.bn);
						cpu.ins_idx = BLD;
					}
					break;
				case 0x0C00:
					if (TEST(ins, 9)) {
						LOG("SBRS r%d, %d", cpu.rd, cpu.bn);
						cpu.ins_idx = SBRS;
					} else {
						LOG("SBRC r%d, %d", cpu.rd, cpu.bn);
						cpu.ins_idx = SBRC;
					}
					break;
				default:
					break;
				}
			}
			if (((ins & 0xF000) == 0x8000) || ((ins & 0xF000) == 0xA000)) {
				if (TEST(ins, 9)) {
					cpu.rd = extract16(ins, 4, 5);
					cpu.imm = ((ins & 0x2000) >> 8) | ((ins & 0x0C00) >> 7) |
								(ins & 0x7);
					if (TEST(ins, 3)) {
						LOG("ST (Y + %d), r%d", cpu.imm, cpu.rd);
						cpu.ins_idx = STDY;
					} else {
						LOG("ST (Z + %d), r%d", cpu.imm, cpu.rd);
						cpu.ins_idx = STDZ;
					}
				} else {
					cpu.rd = extract16(ins, 4, 5);
					cpu.imm = ((ins & 0x2000) >> 8) | ((ins & 0x0C00) >> 7) |
								(ins & 0x7);
					if (TEST(ins, 3)) {
						LOG("LD (Y + %d), r%d", cpu.imm, cpu.rd);
						cpu.ins_idx = LDDY;
					} else {
						LOG("LD (Z + %d), r%d", cpu.imm, cpu.rd);
						cpu.ins_idx = LDDZ;
					}
				}
			}
		}
		break;
	};
}

void decode(uint16_t ins)
{
	unsigned int lz = __builtin_clz(ins) - 16;
	uint8_t rd;
	uint8_t rr;

	if (!ins) {
		LOG("nop");
		return;
	}

	if (ins == 0xFFFF) {
		cpu.ins_idx = BKPT;
		return;
	}

	switch(lz) {
	case 7:
		cpu.rd = extract16(ins, 4, 4);
		cpu.rr = extract16(ins, 0, 4);
		cpu.rd *= 2;
		cpu.rr *= 2;
		cpu.ins_idx = MOVW;
		LOG("MOVW r%d:r%d  r%d:r%d", cpu.rd+1, cpu.rd, cpu.rr+1, cpu.rr);
		break;
	case 6:
		if (TEST(ins, 8) == 0) {
			LOG("muls");
			break;
		}
		if (TEST(ins, 7)) {
			LOG("fmuls");
		} else {
			if (TEST(ins, 3)) {
				LOG("fmul");
			} else {
				LOG("mulsu");
			}
		}
		break;
	case 5:
		cpu.rr = (extract16(ins, 9, 1) << 4 ) | extract16(ins, 0, 4);
		cpu.rd = extract16(ins, 4, 5);
		LOG("CPC r%d r%d", cpu.rd, cpu.rr);
		cpu.ins_idx = CPC;
		break;
	case 4:
		cpu.rr = (extract16(ins, 9, 1) << 4 ) | extract16(ins, 0, 4);
		cpu.rd = extract16(ins, 4, 5);
		if (TEST(ins, 10)) {
			cpu.ins_idx = ADD;
			LOG("ADD r%d r%d", cpu.rd, cpu.rr);
		} else {
			cpu.ins_idx = SBC;
			LOG("SBC r%d r%d", cpu.rd, cpu.rr);
		}
		break;
	case 3:
		cpu.rd = extract16(ins, 4, 5);
		cpu.rr = (extract16(ins, 9, 1) << 4 ) | extract16(ins, 0, 4);
		switch (ins & 0x0C00) {
		case 0x0:
			LOG("cpse r%d r%d", cpu.rd, cpu.rr);
			break;
		case 0x0400:
			LOG("CP r%d r%d", cpu.rd, cpu.rr);
			cpu.ins_idx = CP;
			break;
		case 0x0800:
			LOG("SUB r%d r%d", cpu.rd, cpu.rr);
			cpu.ins_idx = SUB;
			break;
		case 0x0C00:
			LOG("ADC r%d r%d", cpu.rd, cpu.rr);
			cpu.ins_idx = ADC;
			break;
		default:
			break;
		}
		break;
	case 2:
		{
			if (TEST(ins, 12)) {
				cpu.rd = extract16(ins, 4, 4) + 16;
				cpu.imm = ((ins & 0x0F00) >> 4) | (ins & 0xF);
				cpu.ins_idx = CPI;
				LOG("CPI r%d 0x%x", cpu.rd, cpu.imm);
				break;
			}
			cpu.rr = (extract16(ins, 9, 1) << 4 ) | extract16(ins, 0, 4);
			cpu.rd = extract16(ins, 4, 5);
			switch (extract16(ins, 10, 2)) {
			case 0:
				LOG("AND r%d, r%d", cpu.rd, cpu.rr);
				cpu.ins_idx = AND;
				break;
			case 1:
				LOG("EOR r%d, r%d", cpu.rd, cpu.rr);
				cpu.ins_idx = EOR;
				break;
			case 2:
				LOG("OR r%d, r%d", cpu.rd, cpu.rr);
				cpu.ins_idx = OR;
				break;
			case 3:
				LOG("MOV r%d, r%d", cpu.rd, cpu.rr);
				cpu.ins_idx = MOV;
				break;
			default:
				break;
			}
		}
		break;
	case 1:
		switch (ins & 0x3000) {
		case 0x0:
			cpu.rd = extract16(ins, 4, 4) + 16;
			cpu.imm = ((ins & 0x0F00) >> 4) | (ins & 0xF);
			LOG("SBCI r%d, 0x%x", cpu.rd, cpu.imm);
			cpu.ins_idx = SBCI;
			break;
		case 0x1000:
			cpu.rd = extract16(ins, 4, 4) + 16;
			cpu.imm = ((ins & 0x0F00) >> 4) | (ins & 0xF);
			LOG("SUBI r%d, 0x%x", cpu.rd, cpu.imm);
			cpu.ins_idx = SUBI;
			break;
		case 0x2000:
			cpu.rd = extract16(ins, 4, 4) + 16;
			cpu.imm = ((ins & 0x0F00) >> 4) | (ins & 0xF);
			LOG("ORI r%d, 0x%x", cpu.rd, cpu.imm);
			cpu.ins_idx = ORI;
			break;
		case 0x3000:
			cpu.rd = extract16(ins, 4, 4) + 16;
			cpu.imm = ((ins & 0x0F00) >> 4) | (ins & 0xF);
			LOG("ANDI r%d, 0x%x", cpu.rd, cpu.imm);
			cpu.ins_idx = ANDI;
			break;
		default:
			break;
		}
		break;
	default:
		decode_others(ins);
		break;
	};
}

void dump_all() {
	int i;

	printf("===============================================================\n");
	for (i = 0; i < SREG; i++){
		if ((i % 2) == 0)
			printf("\n");
		printf("%s(0x%x): 0x%x  \t", ioreg_names[i], i, cpu.sram.regs.ioreg[i]);
	}
	printf("\nSREG: %c %c %c %c %c %c %c %c",
			TEST(cpu.sram.regs.ioreg[SREG], 7) ? 'I': ' ',
			TEST(cpu.sram.regs.ioreg[SREG], 6) ? 'T': ' ',
			TEST(cpu.sram.regs.ioreg[SREG], 5) ? 'H': ' ',
			TEST(cpu.sram.regs.ioreg[SREG], 4) ? 'S': ' ',
			TEST(cpu.sram.regs.ioreg[SREG], 3) ? 'V': ' ',
			TEST(cpu.sram.regs.ioreg[SREG], 2) ? 'N': ' ',
			TEST(cpu.sram.regs.ioreg[SREG], 1) ? 'Z': ' ',
			TEST(cpu.sram.regs.ioreg[SREG], 0) ? 'C': ' '
			);
	printf("\n");
	for (i = 0; i < 32; i++){
		if ((i % 4) == 0)
			printf("\n");
		printf("r%d: 0x%x \t", i, cpu.sram.regs.reg[i]);
	}
	printf("\n X: 0x%x%x Y: 0x%x%x Z: 0x%x%x\n", cpu.sram.regs.reg[27],
			cpu.sram.regs.reg[26], cpu.sram.regs.reg[29],
			cpu.sram.regs.reg[28], cpu.sram.regs.reg[31], cpu.sram.regs.reg[30]);
	printf("\nPC:0x%x\n", cpu.pc);
	printf("===============================================================\n");
}

static int one_cycle()
{
	/* fetch */
	/* decode */
	/* execute */
	cpu.ins = fetch(cpu.pc);
	cpu.ins_count++;
	decode(cpu.ins);
	cpu.pc += 2;
	return execute();
}

void free_run()
{
	while (1) {
		if (one_cycle() == BKPT)
			break;
	//	delay();
	}
}

static sigset_t set;
static int sig;
static int *sigptr = &sig;
int _COMMAND = -1;
extern uint16_t bp_addr;
extern uint16_t mem_addr;
extern uint16_t nbytes;
extern uint32_t nins;

static void handle_bkpt()
{
	uint16_t ins;
	uint16_t addr = bp_addr & 0xFFF;

	ins = *(uint16_t *)(cpu.rom + addr);
	insert_bkpt(addr, ins);
	*(uint16_t *)(cpu.rom + addr) = 0xFFFF;
}

void under_monitor()
{
	int i;

	while (1) {
		/* wait for signal */
		if (sigwait(&set, sigptr) == -1) {
			perror("sigwait failed\n");
			return;
		} else {
			if(*sigptr != SIGUSR1) {
				printf("sigwait returned with sig: %d\n", *sigptr);
				return;
			}
		}
		/* read cmd from monitor and execute */
		switch (_COMMAND) {
		case SS:
			one_cycle();
			break;
		case NINS:
			while (nins --)
				one_cycle();
			break;
		case CON:
			free_run();
			break;
		case RUN:
			free_run();
			break;
		case BPT:
			handle_bkpt();
			break;
		case DUMP:
			dump_all();
			break;
		case MEM:
			if (nbytes) {
				for (i = 0; i < nbytes; i++)
					printf("%x: 0x%x\n", mem_addr + i,
							cpu.sram.ram[mem_addr + i]);
			} else {
					printf("%x: 0x%x\n", mem_addr, cpu.sram.ram[mem_addr]);
			}
			break;
		default:
			break;
		}
	}
}

void *cpu_fn(void *data)
{
	int i;
	int *halt = data;

	sigemptyset(&set); 
	sigaddset(&set, SIGUSR1);
	sigprocmask(SIG_BLOCK, &set, NULL );

	printf("Atmega8 MCU, created id:%lu\n", pthread_self());
	if (*halt) {
		printf("CPU halted\n");
		cpu.mode = 1;
		log = 1;
		under_monitor();
	} else {
		cpu.mode = 0;
		free_run();
	}
}

void uart_handler(char ch) {
}
