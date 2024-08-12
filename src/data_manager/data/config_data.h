#pragma once
#include <fstream>
#include "object.h"
#include "json/value.h"
#include "json/reader.h"

class JsonParser : public Object
{
public:
    JsonParser(const std::string& file)
    {
        if (!IsJsonFile(file))
        {
            _console->error("{} is not json file!", file);
            return;
        }

        ParseJsonFile(file);
    }

    ~JsonParser() = default;

public:
    template<typename T, typename... Args>
    T Get(std::string key, Args&&... args)
    {
        Json::Value json_node = _json_node;

        auto getNode = [this, &json_node](const auto& arg)
        {
            if (json_node[arg].isObject())
            {
                json_node = json_node[arg];
            }
            else
            {
                _console->error("{} is not a object!", arg);
                return;
            }
        };

        (getNode(std::forward<Args>(args)), ...);

        try
        {
            if (json_node.isMember(key))
            {
                return json_node[key].as<T>();
            }
            else
            {
                throw std::runtime_error("no key named: " + key);
            }
        }
        catch (const std::exception&)
        {
            //log
            _console->error("no key named: {}", key);
            throw;
        }
    }

    auto& GetJsonNode(const std::string& key)
    {
        if (!_json_node[key].isObject())
        {
            _console->warn("Can not find key: {}", key);
        }

        return _json_node[key];
    }

    bool HasNode(const std::string& key)
    {
        return _json_node.isMember(key);
    }

private:

    bool IsJsonFile(const std::string& file)
    {
        auto length = file.length();
        return (length >= 5 && file.substr(length - 5) == ".json");
    }

    void ParseJsonFile(const std::string& file)
    {
        try
        {
            std::ifstream filestream(file);
            if (!filestream.is_open())
            {
                _console->error("failed to open file: {}", file);
                return;
            }

            Json::CharReaderBuilder readerBuilder;
            Json::parseFromStream(readerBuilder, filestream, &_json_node, nullptr);

            //关闭文件
            filestream.close();

            //此处可以删除注释，用于以后做支持注释的扩展
        }
        catch (const std::exception&)
        {
            _console->error("Error parsing JSON");
        }
    }

private:
    Json::Value _json_node;
};