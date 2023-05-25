IDIR  += -I./

BOOST_ROOT = /opt/homebrew/Cellar/boost/1.81.0_1/include
CIMG_ROOT = /opt/homebrew/Cellar/cimg/3.2.4/include
X11_ROOT = /opt/homebrew/Cellar/libx11/1.8.4/include
X11_ROOT_APP = /opt/homebrew/Cellar/xorgproto/2022.2/include/

CC = g++ -I$(BOOST_ROOT) -I$(CIMG_ROOT) -I$(X11_ROOT) -I$(X11_ROOT_APP)

CFLAGS += -std=c++11
CFLAGS += -lfftw3 
CFLAGS += -L/opt/homebrew/Cellar/libx11/1.8.4/lib -lX11
CFLAGS += -L. -lforme

LDFLAGS := -L$(BOOST_ROOT)/stage/lib -lboost_filesystem

program: 
	$(CC) test_fm2.cpp -o main $(CFLAGS) $(IDIR)
	
	
	
