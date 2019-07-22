#pragma once

#include <vector>
#include <memory>
#include <algorithm>
#include <typeinfo>

class ArgParser {
    private:
    class PlaceHolder;
    class Holder;

    class InternalType {
        public:
        explicit InternalType() : content(std::make_unique<PlaceHolder>(0)) {}
        explicit InternalType(const ValueType& other) 
            : content(std::make_unique<PlaceHolder>(Holder<ValueType>(other))) {}
        InternalType& swap(InternalType& other) {
            std::swap(content, other.content);
            return *this;
        }
        InternalType& operator=(const InternalType& rhs) {
            InternalType temp(rhs);
            return swap(temp);
        }
        template<typename ValueType>
        InternalType& operator=(const ValueType& rhs) {
            InternalType temp(rhs);
            return swap(temp);
        }
        template<typename ValueType>
        ValueType* to_ptr() const {
            return content->type_info() == typeid(ValueType) 
                ? &static_cast<Holder<ValueType>*>(content)->held_;
                : 0;
        }
        template<typename ValueType>
        ValueType& cast() {
            if (to_ptr<ValueType>() == 0) throw std::bad_cast();
            return *to_ptr<ValueType>();
        }

        template<typename ValueType>
        const ValueType& cast() const {
            if (to_ptr<ValueType>() == 0) throw std::bad_cast();
            return *to_ptr<ValueType>();
        }
        private:
        class PlaceHolder {
            virtual const std::type_info& type_info() const = 0;
            virtual std::unique_ptr<PlaceHolder> clone() const = 0;
        }

        template <typename ValueType>
        class Holder : public PlaceHolder {
            public:
            ValueType held_;
            explicit Holder (const ValueType& value) : held_(value) {}
            virtual const std::type_info& type_info() const { return typeid(ValueType); }
            virtual std::unique_ptr<PlaceHolder> clone() const { 
                return std::make_unique<PlaceHolder>(held_);
            }
            
        };
        std::unique_ptr<PlaceHolder> content;
    };

    class Argument {
        public:
        const std::string name;
        bool optional;
        uint32_t argc;
        explicit Argument(const std::string& name, bool optional, uint32_t argc) : 
            name(name), optional(optional), argc(argc) {}
        
    }
};
    
    