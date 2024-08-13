#pragma once
#include <fstream>
#include "object.h"
#include "json/value.h"
#include "json/reader.h"

class ConfigData : public Object
{
public:
    /**
     * @brief constructor
     * @param file json config file 
    */
    ConfigData(const std::string& file)
    {
        if (!IsJsonFile(file))
        {
            _console->error("{} is not json file!", file);
            return;
        }

        ParseJsonFile(file);
    }

    ~ConfigData() = default;

public:
    /**
     * @brief get value
     * @tparam T 
     * @tparam ...Args 
     * @param key 
     * @param ...args 
     * @return value
    */
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

    /**
     * @brief get json node
     * @param key 
     * @return json node
    */
    auto& GetJsonNode(const std::string& key)
    {
        if (!_json_node[key].isObject())
        {
            _console->warn("Can not find key: {}", key);
        }

        return _json_node[key];
    }

    /**
     * @brief Check whether there is key Node
     * @param key node name
     * @return true or false
    */
    bool HasNode(const std::string& key)
    {
        return _json_node.isMember(key);
    }

private:

    /**
     * @brief Check whether the file is json file
     * @param file file path
     * @return true or false
    */
    bool IsJsonFile(const std::string& file)
    {
        auto length = file.length();
        return (length >= 5 && file.substr(length - 5) == ".json");
    }

    /**
     * @brief parse json file
     * @param file 
    */
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

            ///close file
            filestream.close();
        }
        catch (const std::exception&)
        {
            _console->error("Error parsing JSON");
        }
    }

private:
    Json::Value _json_node;
};