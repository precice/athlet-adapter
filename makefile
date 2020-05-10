include ../mk.plugin

FC_LIBRARY_DIRS.%    += $(shell pkg-config --variable=libdir libprecice)
FC_LIBRARIES.%       += precice
