/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */

/*! 
 *  \file socket_utils.hpp
 *  \brief Utility Functions for POSIX Sockets
 *  
 *  This file contains a set of utility functions for socket calls, such as
 *  functions to create a connection, to open a new connection, and to
 *  receive and send data through a socket.
 *  
 */

/* ***************************************************************************
 *  
 *  This program is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License version 3 as
 *  published by the Free Software Foundation.
 *  
 *  This program is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 *  License for more details.
 *  
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software Foundation,
 *  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *  
 ****************************************************************************
 */

/*  Authors: Tiziano De Matteis <dematteis@di.unipi.it>, Gabriele Mencagli <mencagli@di.unipi.it>
 *  March 2016
 */

#ifndef _SOCKET_UTILS_H
#define _SOCKET_UTILS_H

// include
#include <stdio.h>
#include <errno.h>
#include <fcntl.h>
#include <netdb.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <arpa/inet.h>

// define
#define BOTTLENECK_ERR -1
#define MAXHOSTNAMELEN 256

using namespace std;

// function to connect_to a given ip address and port
int socket_connect(const char *ip_address, int port) {
	// prepare the sockaddr_in struct
	struct sockaddr_in sockAddress;
	sockAddress.sin_family = AF_INET;
	sockAddress.sin_port = htons(port);
	sockAddress.sin_addr.s_addr = inet_addr(ip_address);
	// if instead of an ip we have a hostname
	struct hostent *hostEntity;
	char hnamebuf[MAXHOSTNAMELEN];
	if(sockAddress.sin_addr.s_addr == (u_int)-1) {
		// try to translate the hostname into an ip address
		hostEntity = gethostbyname(ip_address);
		if(!hostEntity) {
			perror("Error gethostbyname() call");
			exit(-1);
		}
		else {
			sockAddress.sin_family = hostEntity->h_addrtype;
			bcopy(hostEntity->h_addr, (caddr_t)&sockAddress.sin_addr, hostEntity->h_length);
			strncpy(hnamebuf, hostEntity->h_name, sizeof(hnamebuf)-1);
		}
	}
	// opening of the TCP/IP socket
	int s;
	if((s = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0) {
		perror("Error socket() call");
		exit(-1);
	}
	// while loop until a connection has been created
	while(connect(s, (struct sockaddr*) &sockAddress, sizeof(struct sockaddr)) < 0);
	return s;
}

// function to retrieve a number of TCP/IP connections on a given port (returns an array of socket descriptors)
int *socket_accept(int num_conn, int port) {
	if(num_conn > 0) {
		// prepare the sockaddr_in struct
		struct sockaddr_in serverSocketAddress;
		struct sockaddr_in partnerAddr;
		socklen_t sockaddr_len = (socklen_t) sizeof(partnerAddr);
		serverSocketAddress.sin_family = AF_INET;
		serverSocketAddress.sin_port = htons(port);
		serverSocketAddress.sin_addr.s_addr = INADDR_ANY;
		// open server socket
		int serverSocket;
		if((serverSocket = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0) {
			perror("Error socket() call");
			exit(-1);
		}
		// enabling option to reuse the local port without timeout issues
		int yes = 1;
		setsockopt(serverSocket, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(int));
		// bind of the server socket
		if(bind(serverSocket, (struct sockaddr *) &serverSocketAddress, sizeof(serverSocketAddress)) < 0) {
			perror("Error bind() call");
			exit(-1);
		}
		// listen on the server socket
		if(listen(serverSocket, 1) < 0) {
			perror("Error list() call");
			exit(-1);
		}
		// while loop until all the connections have been created
		int *s = (int *) malloc(sizeof(int) * num_conn);
		for(size_t i=0; i<num_conn; i++) {
			if((s[i] = accept(serverSocket, (struct sockaddr *) &partnerAddr, &sockaddr_len)) < 0) {
				perror("Error accept() call");
				exit(-1);
			}
		}
		return s;
	}
	else return NULL;
}

// function to send a message over an open socket
size_t socket_send(int s, void *msg, size_t len) {
	size_t sentBytes = 0;
	size_t sent = 0;
	char *msg_char = (char *) msg;
	// do many sends until all the bytes have been transmitted
	while(sentBytes < len) {
		// sent = send(s, msg_char, len, (int) MSG_DONTWAIT);
		sent = send(s, msg_char, len, 0);
		if(sent == -1 && ((errno == EWOULDBLOCK) || (errno == EAGAIN))) {
			// non blocking mode
			perror("Error send(): buffer full");
			return BOTTLENECK_ERR;
		}
		if(sent == -1) {
			// blocking mode
			perror("Error send() call");
			return 0;
		}
		else if(sent <=  len) {
			sentBytes +=  sent;
			msg_char +=  (sent);
		}
	}
	return sentBytes;
}

// function to receive a message from an open socket
size_t socket_receive(int s, void *vtg, size_t len) {
	char *vtg_char = (char *) vtg;
	// receive
	size_t recv_tot = 0;
	while(recv_tot < len) {
		size_t received = recv(s, vtg_char+recv_tot, len-recv_tot, MSG_WAITALL);
		if(received == -1) {
			perror("Error recv() call");
			return 0;
		}
		if(received == 0) {
			// the client side has gracefully closed the socket
			return 0;
		}
		else recv_tot += received;
	}
	return recv_tot;
}

// non-deterministic selection of more input sockets
int socket_select(fd_set *fds, int *s_array, int nConnections, int lastSocket, int max) {
	fd_set readfds;
	int ready;
	// copy of the socket set
	readfds = *fds;
	// select until we have data in at least one of the socket
	do {
		ready = select(max + 1, &readfds, NULL, NULL, NULL);
	} while(ready == -1);
	// we have found a socket with something
	bool found = false;
	// scan the socket to find the first ready one (fair policy)
	int index = (lastSocket+1) % nConnections;
	int count = 0;
	while((!found) && (count < nConnections)) {
		// check if the socket has data to read
		if(FD_ISSET(s_array[index], &readfds)) found = true;
		else {
			index = (index+1) % nConnections;
			count++;
		}
	}
	if(found) return index;
	else return -1;
}

// close an open socket connection
int socket_close(int socket) {
	return close(socket);
}

// set a socket descriptor in blocking mode or not
bool socket_setblocking(int fd, bool blocking) {
   if (fd < 0) return false;
   else {
   		int flags = fcntl(fd, F_GETFL, 0);
   		if (flags < 0) return false;
   		flags = blocking ? (flags&~O_NONBLOCK):(flags|O_NONBLOCK);
   		return (fcntl(fd, F_SETFL, flags) == 0) ? true:false;
	}
}

#endif
