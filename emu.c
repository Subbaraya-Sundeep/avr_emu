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

#define MAX_CMDS (sizeof(cmds)/sizeof(cmds[0]))

extern int _COMMAND;
uint16_t bp_addr;
uint16_t mem_addr;
uint8_t  nbytes;
uint32_t nins;

extern void *cpu_fn(void *data);
extern uint32_t loadbin(void);

typedef void (*handler)(char ch);

extern void uart_handler(char ch);
void mon_handler(char ch);
int do_help(int argc, char *argv[]);
int do_exit(int argc, char *argv[]);
int do_singlestep(int argc, char *argv[]);
int do_con(int argc, char *argv[]);
int do_run(int argc, char *argv[]);
int do_break(int argc, char *argv[]);
int do_registers(int argc, char *argv[]);
int do_mem(int argc, char *argv[]);
int do_nins(int argc, char *argv[]);

static int focus;
static char line[80]; 
static int ind; 
static pthread_t cpu;

static struct stdin_handler {
	int fd;
	handler rd_handler;
} stdin_handlers[] = {
	{0, uart_handler},
	{0, mon_handler}
};

struct {
	const char *cmd;
	int argc;
	int (*cmd_handler)(int argc, char *argv[]);
	const char *help;
} cmds [] = {
	{"help", 0, do_help, "print usage"},	
	{"ss", 0, do_singlestep, "single stepping"},	
	{"con", 0, do_con, "continue execution"},	
	{"run", 0, do_run, "run program"},	
	{"break", 0, do_break, "break point"},	
	{"registers", 0, do_registers, "dump registers"},	
	{"x", 0, do_mem, "print content at memory address"},	
	{"nins", 0, do_nins, "run till n instructions"},	
	{"quit", 0, do_exit, "quit program"}
};

int do_help(int argc, char *argv[])
{
	int i;

	printf("\n== Atmega8 emulator and debugger ==\n");
	for (i = 0 ; i < MAX_CMDS; i++)
		printf("\n%s - %s\n", cmds[i].cmd, cmds[i].help);
}

int do_exit(int argc, char *argv[]) {
	exit(0);
}

int do_singlestep(int argc, char *argv[]) {
	_COMMAND = SS;
	pthread_kill(cpu, SIGUSR1);
}

int do_run(int argc, char *argv[]) {
	_COMMAND = RUN;
	pthread_kill(cpu, SIGUSR1);
}

int do_con(int argc, char *argv[]) {
	_COMMAND = CON;
	pthread_kill(cpu, SIGUSR1);
}

int do_break(int argc, char *argv[]) {
	bp_addr = strtoul(argv[1], NULL, 16);
	_COMMAND = BPT;
	pthread_kill(cpu, SIGUSR1);
}

int do_registers(int argc, char *argv[]) {
	_COMMAND = DUMP;
	pthread_kill(cpu, SIGUSR1);
}

int do_mem(int argc, char *argv[]) {
	mem_addr = strtoul(argv[1], NULL, 16);
	if (argc > 1)
		nbytes = strtoul(argv[2], NULL, 16);
	else
		nbytes = 1;
	_COMMAND = MEM;
	pthread_kill(cpu, SIGUSR1);
}

int do_nins(int argc, char *argv[]) {
	nins = strtoul(argv[1], NULL, 10);
	_COMMAND = NINS;
	pthread_kill(cpu, SIGUSR1);
}

static void muxer(int ch)
{
	switch (ch) {
	case 0xD: /* ~ */
	   /* switch b/w MCU uart and monitor */
		break;
	case 0x7e: /* ~ */
	   /* switch b/w MCU uart and monitor */
		focus = !focus;
		/* fall through */
	default:
		stdin_handlers[focus].rd_handler(ch);
		break;
	}
}

void process(char *cmdp, int len)
{
	char delim[2] = " ";
	char *token;
	int argc = 0, i;
	char *argv[10];

	if (len == 0)
		return;
	/* get the first token */
	argv[argc] = token = strtok(cmdp, delim);

	/* walk through other tokens */
	while(token != NULL) {
		token = strtok(NULL, delim);
		if (token != NULL)
			argv[++argc] = token;
	}

	for (i = 0; i < MAX_CMDS; i++)
	{
		if (strcmp(argv[0], cmds[i].cmd) == 0)
			cmds[i].cmd_handler(argc, argv);
	}
}

void mon_handler(char ch) 
{
	int len;

	if (!ind) {
		printf("(monitor) ");
		fflush(stdout);
	}
	if (ch == '\n') {
		line[ind] = '\0';
		len = ind;
		ind = 0;
		goto process_cmd;
	} else {
		line[ind++] = ch;
		return;
	}
process_cmd:
	process(line, len);
}

int main(int argc, char* argv[]) {
	fd_set rfds;
	int retval;
	char ch;
	int halt = 0;
	int *haltp = &halt;

	if ((argc > 1) && (strcmp(argv[1], "-D") == 0))
		halt = 1;

	loadbin();
	/* Watch stdin (fd 0) to see when it has input. */
	FD_ZERO(&rfds);
	FD_SET(0, &rfds);
	memset(line, 0, sizeof(line));

	if(pthread_create(&cpu, NULL, cpu_fn, (void *)haltp)) {
		fprintf(stderr, "Error creating thread\n");
		return -1;
	}

	while (1) {
		retval = select(1, &rfds, NULL, NULL, NULL);
		if (retval == -1) {
			perror("select()");
			return -1;
		}

		if (FD_ISSET(0, &rfds)) {
			read(0, &ch, 1);
			muxer(ch);
		}
	}

	if(pthread_join(cpu, NULL)) {
		fprintf(stderr, "Error joining thread\n");
	}

	return 0;
}
