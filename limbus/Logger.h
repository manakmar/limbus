#pragma once

#ifndef LIMBUS_LOGGER_H
#define LIMBUS_LOGGER_H

#include <limbus/Config.h>

#define LOGE(msg) std::cerr << msg << std::endl

#ifdef LOG_ENABLED
	#include <iostream>
	#define LOG(msg) std::cout << msg << std::endl

#ifdef _DEBUG
	#define LOGD(msg) std::cout << msg << std::endl
#else
	#define LOGD(msg) ()
#endif

#else
	#define LOG(msg) 
	#define LOGD(msg)
#endif

#endif
