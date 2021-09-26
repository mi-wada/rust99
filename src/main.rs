fn main() {}

#[cfg(test)]
mod tests {
    use std::{collections::HashSet, hash::Hash};

    #[test]
    fn p01() {
        let li = vec![1, 1, 2, 3, 5, 8];
        fn last(v: &Vec<i32>) -> Option<i32> {
            v.last().copied()
        }
        assert_eq!(Some(8), last(&li));
    }

    #[test]
    fn p02() {
        let li = vec![1, 1, 2, 3, 5, 8];
        fn penultimate(v: &Vec<i32>) -> Option<i32> {
            v[..v.len() - 1].last().copied()
        }
        assert_eq!(Some(5), penultimate(&li));
    }

    #[test]
    fn p03() {
        let li = vec![1, 1, 2, 3, 5, 8];
        fn nth(v: &Vec<i32>, n: usize) -> Option<i32> {
            v.iter().nth(n).copied()
        }
        assert_eq!(Some(2), nth(&li, 2));
    }

    #[test]
    fn p04() {
        let li = vec![1, 1, 2, 3, 5, 8];
        fn length(v: &Vec<i32>) -> usize {
            v.len()
        }
        assert_eq!(6, length(&li));
    }

    #[test]
    fn p05() {
        let li = vec![1, 1, 2, 3, 5, 8];
        fn reverse(v: &Vec<i32>) -> Vec<i32> {
            let mut res = Vec::with_capacity(v.len());
            let mut copied_v = v.clone();
            while let Some(x) = copied_v.pop() {
                res.push(x);
            }
            res
        }
        assert_eq!(vec![8, 5, 3, 2, 1, 1], reverse(&li));
    }

    #[test]
    fn p06() {
        fn is_palindrome(v: &Vec<i32>) -> bool {
            for i in 0..v.len() / 2 {
                if v[i] != v[v.len() - 1 - i] {
                    return false;
                }
            }
            true
        }
        let palindrome = vec![1, 2, 3, 2, 1];
        assert_eq!(true, is_palindrome(&palindrome));
        let not_palindrome = vec![1, 3, 3, 2, 1];
        assert_eq!(false, is_palindrome(&not_palindrome));
    }

    #[test]
    fn p07() {
        #[derive(PartialEq, Debug)]
        struct List<T: Copy>(Vec<ListItem<T>>);

        impl<T: Copy> List<T> {
            fn new() -> List<T> {
                List(vec![])
            }

            fn from(items: Vec<ListItem<T>>) -> List<T> {
                List(items)
            }

            fn push(&mut self, element: ListItem<T>) {
                self.0.push(element)
            }

            fn append(&mut self, other: List<T>) {
                for x in other.0.into_iter() {
                    self.0.push(x);
                }
            }
        }

        #[derive(Debug, PartialEq)]
        enum ListItem<T: Copy> {
            Item(T),
            NestedItem(Box<List<T>>),
        }

        fn item<T: Copy>(x: T) -> ListItem<T> {
            ListItem::Item(x)
        }

        fn nested_item<T: Copy>(li: List<T>) -> ListItem<T> {
            ListItem::NestedItem(Box::new(li))
        }

        let not_flat_list = List::from(vec![
            item(10),
            nested_item(List::from(vec![item(20), item(30)])),
        ]);

        fn flatten<T: Copy>(li: List<T>) -> List<T> {
            let mut res = List::new();
            for x in li.0 {
                match x {
                    ListItem::Item(x) => res.push(item(x)),
                    ListItem::NestedItem(li) => res.append(flatten(*li)),
                }
            }
            res
        }

        assert_eq!(
            flatten(not_flat_list),
            List::from(vec![item(10), item(20), item(30)])
        );
    }

    #[test]
    fn p08() {
        let li = vec![
            'a', 'a', 'a', 'a', 'b', 'c', 'c', 'a', 'a', 'd', 'e', 'e', 'e', 'e',
        ];
        fn compress<T>(v: &Vec<T>) -> Vec<T>
        where
            T: Copy + Eq + Hash,
        {
            let mut repeated_detecter = HashSet::<T>::new();
            let mut compressd = vec![];
            for x in v {
                if !repeated_detecter.contains(x) {
                    compressd.push(*x);
                }
                repeated_detecter.clear();
                repeated_detecter.insert(*x);
            }
            compressd
        }
        assert_eq!(compress(&li), vec!['a', 'b', 'c', 'a', 'd', 'e']);
    }

    #[test]
    fn p09() {
        let li = vec![
            'a', 'a', 'a', 'a', 'b', 'c', 'c', 'a', 'a', 'd', 'e', 'e', 'e', 'e',
        ];
        fn pack<T: Copy + PartialEq>(v: &Vec<T>) -> Vec<Vec<T>> {
            let mut packed = vec![vec![]];
            if v.is_empty() {
                return packed;
            }
            let mut prev_value = &v[0];
            for x in v {
                if x == prev_value {
                    packed.last_mut().unwrap().push(*x);
                } else {
                    packed.push(vec![*x]);
                }
                prev_value = x;
            }
            packed
        }
        assert_eq!(
            pack(&li),
            vec![
                vec!['a', 'a', 'a', 'a'],
                vec!['b'],
                vec!['c', 'c'],
                vec!['a', 'a'],
                vec!['d'],
                vec!['e', 'e', 'e', 'e']
            ]
        );
    }

    #[test]
    fn p10() {
        let li = vec![
            'a', 'a', 'a', 'a', 'b', 'c', 'c', 'a', 'a', 'd', 'e', 'e', 'e', 'e',
        ];
        fn pack<T: Copy + PartialEq>(v: &Vec<T>) -> Vec<Vec<T>> {
            let mut packed = vec![vec![]];
            if v.is_empty() {
                return packed;
            }
            let mut prev_value = &v[0];
            for x in v {
                if x == prev_value {
                    packed.last_mut().unwrap().push(*x);
                } else {
                    packed.push(vec![*x]);
                }
                prev_value = x;
            }
            packed
        }
        fn encode<T: Copy + PartialEq>(v: &Vec<T>) -> Vec<(usize, T)> {
            let mut encoded = vec![];
            for x in pack(v) {
                encoded.push((x.len(), *x.last().unwrap()));
            }
            encoded
        }
        assert_eq!(
            encode(&li),
            vec![(4, 'a'), (1, 'b'), (2, 'c'), (2, 'a'), (1, 'd'), (4, 'e')]
        );
    }

    fn p12() {
        let li = vec![(4, 'a'), (1, 'b'), (2, 'c'), (2, 'a'), (1, 'd'), (4, 'e')];
        fn decode<T: Copy + PartialEq>(v: &Vec<(usize, T)>) -> Vec<T> {
            let mut decoded = vec![];
            for x in v {
                for _ in 0..x.0 {
                    decoded.push(x.1);
                }
            }
            decoded
        }
        assert_eq!(
            decode(&li),
            vec!['a', 'a', 'a', 'a', 'b', 'c', 'c', 'a', 'a', 'd', 'e', 'e', 'e', 'e']
        );
    }

    #[test]
    fn p13() {
        let li = vec![
            'a', 'a', 'a', 'a', 'b', 'c', 'c', 'a', 'a', 'd', 'e', 'e', 'e', 'e',
        ];
        fn encode<T: Copy + PartialEq>(v: &Vec<T>) -> Vec<(usize, T)> {
            if v.is_empty() {
                vec![]
            } else {
                let mut encoded = vec![(1, v[0].clone())];
                for x in v[1..].iter() {
                    if x == &encoded.last().unwrap().1 {
                        encoded.last_mut().unwrap().0 += 1;
                    } else {
                        encoded.push((1, *x));
                    }
                }
                encoded
            }
        }
        assert_eq!(
            encode(&li),
            [(4, 'a'), (1, 'b'), (2, 'c'), (2, 'a'), (1, 'd'), (4, 'e')]
        );
    }
}
