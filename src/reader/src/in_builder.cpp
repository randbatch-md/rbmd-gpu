#include "../include/in_builder.h"
#include "model/md_data.h"

InBuilder::InBuilder():
	_md_data(DataManager::getInstance().getMDData())
{}