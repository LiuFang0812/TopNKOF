package util;

public class LinkQueue<E> {
    // 链栈的节点
    public class ListNode<E> {
        E e;
        ListNode<E> next;

        public ListNode() {
        }

        public ListNode(E e, ListNode next) {
            this.e = e;
            this.next = next;
        }

		public E getE() {
			return e;
		}

		public void setE(E e) {
			this.e = e;
		}

		public ListNode<E> getNext() {
			return next;
		}

		public void setNext(ListNode<E> next) {
			this.next = next;
		}
        
        
    }
    
    private ListNode front;// 队列头，允许删除  
    private ListNode rear;// 队列尾，允许插入  
    private int size; //队列当前长度 
    
    public LinkQueue() {
        front = null;
        rear = null;
    }
    
    //判空
      public boolean empty(){
          return size==0;
      }
      
      //插入
      public boolean add(E e){
          if(empty()){    //如果队列为空
              front = new ListNode(e,null);//只有一个节点，front、rear都指向该节点
              rear = front;
          }else{
        	  ListNode<E> newNode = new ListNode<E>(e, null);
              rear.next = newNode; //让尾节点的next指向新增的节点
              rear = newNode; //以新节点作为新的尾节点
          }
          size ++;
          return true;
      }
      
      //返回队首元素，但不删除
      public ListNode<E> peek(){
          if(empty()){
              throw new RuntimeException("空队列异常！");
          }else{
              return front;
          }
      }
      
      //出队
      public ListNode<E> poll(){
          if(empty()){
              throw new RuntimeException("空队列异常！");
          }else{
        	  ListNode<E> value = front; //得到队列头元素
              front = front.next;//让front引用指向原队列头元素的下一个元素
              value.next = null; //释放原队列头元素的next引用
              size --;
              return value;
          }        
      }
      
      //队列长度
      public int length(){
          return size;
      }

	public ListNode getFront() {
		return front;
	}

	public void setFront(ListNode front) {
		this.front = front;
	}

	public ListNode getRear() {
		return rear;
	}

	public void setRear(ListNode rear) {
		this.rear = rear;
	}


	
      
}
